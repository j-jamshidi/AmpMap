# app.py
from flask import Flask, render_template, jsonify, send_file, abort
import os
import glob
import re

app = Flask(__name__)

# Base directory for all operations
BASE_DIR = '/Users/javadjamshidi/Desktop/Runs'

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
        
        # Look through barcode01 to barcode24
        for i in range(1, 25):
            barcode_dir = os.path.join(result_dir, f'barcode{i:02d}')
            if not os.path.exists(barcode_dir):
                continue
                
            # Find all report files
            report_files = glob.glob(os.path.join(barcode_dir, '*_report.txt'))
            for report_file in report_files:
                # Extract episode name from report filename
                episode = os.path.basename(report_file).replace('_report.txt', '')
                if episode:  # Only add if episode name is not empty
                    samples.add(episode)
        
        return jsonify(sorted(list(samples)))
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/runs/<run_id>/samples/<sample_id>/report')
def get_sample_report(run_id, sample_id):
    """Get report content for a specific sample"""
    try:
        result_dir = os.path.join(BASE_DIR, run_id, 'result')
        report_content = None
        
        # Search through barcode01 to barcode24
        for i in range(1, 25):
            barcode_dir = os.path.join(result_dir, f'barcode{i:02d}')
            if not os.path.exists(barcode_dir):
                continue
                
            report_path = os.path.join(barcode_dir, f'{sample_id}_report.txt')
            if os.path.exists(report_path):
                with open(report_path, 'r') as f:
                    report_content = f.read()
                break
        
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
        
        # Define file name based on type
        file_names = {
            'xml': f'{sample_id}.xml',
            'HapCUT2': f'HapCUT2.log',
            'WhatsHap': f'whatshap.log'
        }
        
        if file_type not in file_names:
            return jsonify({'error': 'Invalid file type'}), 400
            
        file_name = file_names[file_type]
        
        # Search through barcode01 to barcode24
        for i in range(1, 25):
            barcode_dir = os.path.join(result_dir, f'barcode{i:02d}')
            if not os.path.exists(barcode_dir):
                continue
                
            temp_path = os.path.join(barcode_dir, file_name)
            if os.path.exists(temp_path):
                file_path = temp_path
                break
        
        if file_path is None:
            return jsonify({'error': 'File not found'}), 404
            
        return send_file(file_path, as_attachment=True)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run(debug=True, port=5000)