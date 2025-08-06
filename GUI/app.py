from flask import Flask, render_template, jsonify, send_file, abort
import os
import glob
import re
import sqlite3

app = Flask(__name__)

# Base directory for all operations
BASE_DIR = '/Users/javadjamshidi/Desktop/Runs'
DB_PATH = '/Users/javadjamshidi/Desktop/ONT_amplicon_phase/GUI/amplicon_data.db'

def init_db():
    """Initialize the database."""
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS amplicon_data (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            run_id TEXT,
            sample_id TEXT,
            amplicon_size INTEGER,
            total_reads INTEGER,
            passing_qc_reads INTEGER
        )
    ''')
    conn.commit()
    conn.close()

def insert_amplicon_data(run_id, sample_id, amplicon_size, total_reads, passing_qc_reads):
    """Insert amplicon data into the database."""
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()
    cursor.execute('''
        INSERT INTO amplicon_data (run_id, sample_id, amplicon_size, total_reads, passing_qc_reads)
        VALUES (?, ?, ?, ?, ?)
    ''', (run_id, sample_id, amplicon_size, total_reads, passing_qc_reads))
    conn.commit()
    conn.close()

def parse_report_file(report_file, run_id, sample_id):
    """Parse the report file and extract amplicon data."""
    amplicon_size = 0
    total_reads = 0
    passing_qc_reads = 0

    with open(report_file, 'r') as f:
        for line in f:
            if line.startswith("Amplicon length:"):
                amplicon_size = int(line.split(":")[1].strip().split()[0].replace(',', ''))
            elif line.startswith("Total reads:"):
                total_reads = int(line.split(":")[1].strip().split()[0].replace(',', ''))
            elif line.startswith("Passing QC reads:"):
                passing_qc_reads = int(line.split(":")[1].strip().split()[0].replace(',', ''))

    insert_amplicon_data(run_id, sample_id, amplicon_size, total_reads, passing_qc_reads)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/api/runs')
def get_runs():
    """Get list of all runs"""
    try:
        runs = [d for d in os.listdir(BASE_DIR) 
                if os.path.isdir(os.path.join(BASE_DIR, d))]
        return jsonify(sorted(runs))
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/runs/<run_id>/html')
def get_run_html(run_id):
    """Get HTML report for a specific run"""
    try:
        run_dir = os.path.join(BASE_DIR, run_id)
        html_files = glob.glob(os.path.join(run_dir, '*.html'))
        if not html_files:
            return jsonify({'error': 'HTML file not found'}), 404
        
        with open(html_files[0], 'r') as f:
            content = f.read()
        return content
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/runs/<run_id>/samples')
def get_samples(run_id):
    """Get list of samples for a specific run"""
    try:
        samples = set()  # Use set to avoid duplicates
        result_dir = os.path.join(BASE_DIR, run_id, 'result')
        
        # Look through barcode01 to 92
        for i in range(1, 93):
            barcode_dir = os.path.join(result_dir, f'barcode{i:02d}')
            if not os.path.exists(barcode_dir):
                continue
                
            # Find all report files
            report_files = glob.glob(os.path.join(barcode_dir, '*_report.txt'))
            for report_file in report_files:
                # Extract episode name from report filename and add barcode suffix
                episode = os.path.basename(report_file).replace('_report.txt', '')
                if episode:  # Only add if episode name is not empty
                    samples.add(f"{episode}_b{i:02d}")
        
        return jsonify(sorted(list(samples)))
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/runs/<run_id>/samples/<sample_id>/report')
def get_sample_report(run_id, sample_id):
    """Get report content for a specific sample"""
    try:
        result_dir = os.path.join(BASE_DIR, run_id, 'result')
        report_content = None
        
        # Extract the original sample name and barcode number
        match = re.match(r"(.+)_b(\d{2})$", sample_id)
        if not match:
            return jsonify({'error': 'Invalid sample ID format'}), 400
        
        original_sample_id, barcode = match.groups()
        
        # Search through the specified barcode directory
        barcode_dir = os.path.join(result_dir, f'barcode{barcode}')
        report_path = os.path.join(barcode_dir, f'{original_sample_id}_report.txt')
        if os.path.exists(report_path):
            with open(report_path, 'r') as f:
                report_content = f.read()
        
        if report_content is None:
            return jsonify({'error': 'Report not found'}), 404
            
        return report_content
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/runs/<run_id>/samples/<sample_id>/download/<file_type>')
def download_file(run_id, sample_id, file_type):
    """Download specific file types"""
    try:
        result_dir = os.path.join(BASE_DIR, run_id, 'result')
        file_path = None
        
        # Extract the original sample name and barcode number
        match = re.match(r"(.+)_b(\d{2})$", sample_id)
        if not match:
            return jsonify({'error': 'Invalid sample ID format'}), 400
        
        original_sample_id, barcode = match.groups()
        
        # Define file name based on type
        file_names = {
            'xml': f'{original_sample_id}.xml',
            'HapCUT2': f'HapCUT2.log',
            'WhatsHap': f'whatshap.log'
        }
        
        if file_type not in file_names:
            return jsonify({'error': 'Invalid file type'}), 400
            
        file_name = file_names[file_type]
        
        # Search through the specified barcode directory
        barcode_dir = os.path.join(result_dir, f'barcode{barcode}')
        temp_path = os.path.join(barcode_dir, file_name)
        if os.path.exists(temp_path):
            file_path = temp_path
        
        if file_path is None:
            return jsonify({'error': 'File not found'}), 404
            
        return send_file(file_path, as_attachment=True)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/amplicon_data')
def get_amplicon_data():
    """Get amplicon data for visualization"""
    try:
        conn = sqlite3.connect(DB_PATH)
        cursor = conn.cursor()
        cursor.execute('SELECT * FROM amplicon_data')
        data = cursor.fetchall()
        conn.close()
        return jsonify(data)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    init_db()
    app.run(debug=True, host='0.0.0.0', port=5001)