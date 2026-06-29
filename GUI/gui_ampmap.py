from flask import Flask, render_template, jsonify, send_file, abort
import os
import glob
import re
import paramiko
import shlex
from datetime import datetime
from html import escape
from io import BytesIO
from pathlib import Path
from time import sleep
from os.path import basename, relpath
from threading import Thread, Lock
import logging
from reportlab.lib import colors
from reportlab.lib.enums import TA_RIGHT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas as pdfcanvas
from reportlab.platypus import Paragraph, SimpleDocTemplate, Spacer

# Import configuration
from config import HOSTNAME, USERNAME, BASE_PATH, LOCAL_PATH, PEM_PATH, FLASK_HOST, FLASK_PORT

class RemoteFileMonitor:
    def __init__(self, hostname=HOSTNAME, base_path=BASE_PATH, local_path=LOCAL_PATH, pem_path=PEM_PATH):
        self.hostname = hostname
        self.base_path = base_path
        self.local_path = Path(local_path)
        self.pem_path = Path(pem_path)
        self.host_ip = hostname
        self.username = "ubuntu"
        
        # Set up logging
        self.logger = logging.getLogger('RemoteMonitor')
        self.logger.setLevel(logging.INFO)
        
        if not self.pem_path.exists():
            raise FileNotFoundError(f"PEM file not found: {self.pem_path}")
            
        # Create local directory if it doesn't exist
        self.local_path.mkdir(parents=True, exist_ok=True)
        
        # Set up SSH client with proper host key policy
        self.ssh = paramiko.SSHClient()
        self.ssh.load_system_host_keys()
        if '3.24.162.88' in hostname:
            self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        
        self.sftp = None
        self.processed_signals = set()
        self._lock = Lock()

    def connect(self):
        """Establish SSH connection"""
        try:
            if not os.path.exists(self.pem_path):
                raise FileNotFoundError(f"PEM file not found at {self.pem_path}")
                
            # Close existing connections if any
            if self.sftp:
                try:
                    self.sftp.close()
                except:
                    pass
            if self.ssh:
                try:
                    self.ssh.close()
                except:
                    pass
                    
            self.ssh.connect(
                hostname=self.host_ip,
                username=USERNAME,
                key_filename=str(self.pem_path),
                timeout=10
            )
            self.sftp = self.ssh.open_sftp()
            self.logger.info(f"Connected to {self._sanitize_log_input(self.hostname)}")
        except Exception as e:
            self.logger.error(f"Connection failed: {self._sanitize_log_input(str(e))}")
            self.sftp = None
            raise

    def _sanitize_path(self, path):
        """Sanitize path to prevent traversal attacks"""
        return re.sub(r'\.\.[\/\\]', '', path)
    
    def _sanitize_log_input(self, text):
        """Sanitize input for logging to prevent log injection"""
        if isinstance(text, str):
            return text.replace('\n', ' ').replace('\r', ' ')
        return str(text)

    def check_new_signal_files(self):
        """Check for new .analysed signal files"""
        try:
            cmd = f'find {shlex.quote(self.base_path)} -maxdepth 1 -type f -name "*.analysed"'
            _, stdout, stderr = self.ssh.exec_command(cmd)
            
            new_signals = set(stdout.read().decode().splitlines())
            signals_to_process = new_signals - self.processed_signals
            
            if not signals_to_process:
                self.logger.info("No new .analysed files to process")
                return
            
            for signal_file in signals_to_process:
                folder_name = basename(signal_file).replace('.analysed', '')
                folder_name = re.sub(r'[^\w\-_.]', '', folder_name)
                self.process_data_folder(folder_name)
                with self._lock:
                    self.processed_signals.add(signal_file)
                
        except Exception as e:
            self.logger.error(f"Error checking for new signal files: {self._sanitize_log_input(str(e))}")
            raise

    def process_data_folder(self, folder_name):
        """Process a data folder when its signal file is detected"""
        try:
            safe_folder_name = shlex.quote(folder_name)
            data_folder = f"{self.base_path}/{safe_folder_name}"
            
            # Check if the data folder exists
            test_cmd = f'test -d {shlex.quote(data_folder)} && echo "exists"'
            _, stdout, stderr = self.ssh.exec_command(test_cmd)
            if not stdout.read().decode().strip():
                self.logger.error(f"Data folder {self._sanitize_log_input(data_folder)} does not exist")
                return
            
            # Download .html, .csv files from the main folder
            self.download_files_by_pattern(data_folder, ['*.html', '*.csv'])
            
            # Download .log files from result folder
            result_folder = f"{data_folder}/result"
            self.download_files_by_pattern(result_folder, ['*.log'])
            
            # Download .log files from barcode directories
            for barcode in range(1, 93):
                barcode_dir = f"{data_folder}/result/barcode{barcode:02d}"
                self.download_files_by_pattern(barcode_dir, ['*.log', '*.txt', '*.xml'])
            
            # Rename .analysed file to .done
            analysed_file = f"{self.base_path}/{safe_folder_name}.analysed"
            done_file = f"{self.base_path}/{safe_folder_name}.done"
            mv_cmd = f'mv {shlex.quote(analysed_file)} {shlex.quote(done_file)}'
            self.ssh.exec_command(mv_cmd)
            
        except Exception as e:
            self.logger.error(f"Error processing folder {self._sanitize_log_input(folder_name)}: {self._sanitize_log_input(str(e))}")
            raise

    def download_files_by_pattern(self, remote_dir, patterns):
        """Download files matching specific patterns from a remote directory"""
        try:
            _, stdout, stderr = self.ssh.exec_command(f'test -d {shlex.quote(remote_dir)} && echo "exists"')
            if not stdout.read().decode().strip():
                self.logger.info(f"Directory {self._sanitize_log_input(remote_dir)} does not exist, skipping")
                return
            
            safe_patterns = [shlex.quote(p) for p in patterns if re.match(r'^[\w*.-]+$', p)]
            if not safe_patterns:
                return
                
            pattern_str = ' -o '.join([f'-name {p}' for p in safe_patterns])
            cmd = f'find {shlex.quote(remote_dir)} -maxdepth 1 -type f \\( {pattern_str} \\)'
            _, stdout, stderr = self.ssh.exec_command(cmd)
            
            files = stdout.read().decode().splitlines()
            for remote_file in files:
                try:
                    if '..' in remote_file or not remote_file.startswith(self.base_path):
                        self.logger.warning(f"Skipping suspicious file path: {self._sanitize_log_input(remote_file)}")
                        continue
                        
                    relative_path = relpath(remote_file, self.base_path)
                    local_file = self.local_path / relative_path
                    local_file.parent.mkdir(parents=True, exist_ok=True)
                    
                    if self.sftp is None:
                        self.sftp = self.ssh.open_sftp()
                    self.sftp.get(remote_file, str(local_file))
                    self.logger.info(f"Successfully downloaded {self._sanitize_log_input(remote_file)}")
                    
                except Exception as e:
                    self.logger.error(f"Error downloading {self._sanitize_log_input(remote_file)}: {self._sanitize_log_input(str(e))}")
                    # Try to reconnect if SFTP is failing
                    try:
                        self.sftp = self.ssh.open_sftp()
                    except:
                        pass
                    continue
                    
        except Exception as e:
            self.logger.error(f"Error processing directory {self._sanitize_log_input(remote_dir)}: {self._sanitize_log_input(str(e))}")
            raise

    def run(self, check_interval=300):

        self.logger.info(f"Starting monitoring of {self._sanitize_log_input(self.hostname)}:{self._sanitize_log_input(self.base_path)}")
        
        connected = False
        retry_count = 0
        
        while True:
            try:
                if not connected:
                    self.connect()
                    connected = True
                    retry_count = 0
                    
                self.check_new_signal_files()
                sleep(check_interval)
                
            except KeyboardInterrupt:
                self.logger.info("Monitoring stopped by user")
                break
            except Exception as e:
                self.logger.error(f"Unexpected error: {self._sanitize_log_input(str(e))}")
                connected = False
                retry_count += 1
                
                if self.ssh:
                    try:
                        self.ssh.close()
                    except:
                        pass
                        
                if self.sftp:
                    try:
                        self.sftp.close()
                    except:
                        pass
                        
                wait_time = min(check_interval * (2 ** min(retry_count, 3)), 600)
                sleep(wait_time)

class GUI_AmpMap:
    def __init__(self):
        # Set up logging first
        self.logger = logging.getLogger('GUI_AmpMap')
        
        self.logger.info("Initializing GUI AmpMap")
        self.logger.info(f"Local path: {LOCAL_PATH}")
        self.logger.info(f"Remote host: {HOSTNAME}")
        
        # Initialize Flask app
        self.app = Flask(__name__)
        self.init_routes()
        self.logger.info("Flask app initialized")
        
        # Initialize remote monitor
        self.remote_monitor = RemoteFileMonitor()
        self.logger.info("Remote monitor initialized")

    def _get_ampmap_version(self, result_dir):
        """Extract AmpMap version from the run log."""
        ampmap_version = "AmpMap"
        try:
            log_files = glob.glob(os.path.join(result_dir, '*.log'))
            if log_files:
                with open(log_files[0], 'r') as log_file:
                    for line in log_file.readlines()[:5]:
                        if 'AmpMap v' in line:
                            version_match = re.search(r'AmpMap v[\d.]+', line)
                            if version_match:
                                ampmap_version = version_match.group()
                            break
        except Exception as e:
            self.logger.warning(f"Version extraction failed: {str(e)}")
        return ampmap_version

    def _get_sample_report_data(self, run_id, sample_id):
        """Load a sample report and return the data needed by UI/download routes."""
        result_dir = os.path.join(str(LOCAL_PATH), run_id, 'result')

        match = re.match(r"(.+)_b(\d{2})$", sample_id)
        if not match:
            return {'error': 'Invalid sample ID format'}, 400

        original_sample_id, barcode = match.groups()
        barcode_dir = os.path.join(result_dir, f'barcode{barcode}')
        report_path = os.path.join(barcode_dir, f'{original_sample_id}_report.txt')

        if not os.path.exists(report_path):
            return {'error': 'Report not found'}, 404

        with open(report_path, 'r') as f:
            report_content = f.read()

        ampmap_version = self._get_ampmap_version(result_dir)
        report_lines = report_content.split('\n')
        displayed_sample_name = original_sample_id

        for i, line in enumerate(report_lines):
            if 'Report for' in line:
                sample_match = re.search(r'Report for[:\s]*([^=]+)', line)
                if sample_match:
                    displayed_sample_name = sample_match.group(1).strip()
                    report_lines[i] = (
                        f"Report for: {displayed_sample_name}"
                        f"<span class='version'>Pipeline Version: {ampmap_version}</span>"
                    )
                break

        mtime = os.path.getmtime(report_path)
        analysis_datetime = datetime.fromtimestamp(mtime).strftime('%d %b %Y')

        return {
            'result_dir': result_dir,
            'original_sample_id': original_sample_id,
            'displayed_sample_name': displayed_sample_name,
            'barcode': barcode,
            'ampmap_version': ampmap_version,
            'analysis_datetime': analysis_datetime,
            'report_content': '\n'.join(report_lines)
        }, 200

    def _normalize_report_title(self, title):
        return re.sub(r'[^a-z0-9]+', ' ', re.sub(r'=+', '', title or '').lower()).strip()

    def _split_sample_report_sections(self, report_content):
        sections = []
        current_section = []

        for line in report_content.split('\n'):
            line = line.strip()
            if line.startswith('='):
                if current_section:
                    sections.append(current_section)
                current_section = [line]
            else:
                current_section.append(line)

        if current_section:
            sections.append(current_section)

        return sections

    def _build_sample_report_section_objects(self, report_content):
        """Mirror the browser's Sample Report tab section reshaping."""
        sections = self._split_sample_report_sections(report_content)
        first_section = sections[0] if sections else []
        first_content_lines = [line for line in first_section[1:] if line.strip() and '---' not in line]

        section_objs = []
        for section in sections[1:]:
            header = section[0] if section else ''
            title = re.sub(r'=+', '', header).strip()
            section_objs.append({
                'header': header,
                'title': title,
                'titleNorm': self._normalize_report_title(title),
                'lines': section[1:]
            })

        idx_result = next((i for i, s in enumerate(section_objs) if s['titleNorm'] in ('result', 'results')), -1)
        idx_qc_for_move = next((i for i, s in enumerate(section_objs) if s['titleNorm'] == 'quality control'), -1)

        if idx_result != -1 and idx_qc_for_move != -1:
            result_lines = section_objs[idx_result]['lines']
            qc_lines = section_objs[idx_qc_for_move]['lines']

            detailed_cat_start = next((i for i, line in enumerate(result_lines) if 'Detailed categorisation of reads:' in line), -1)
            if detailed_cat_start != -1:
                detailed_cat_end = next(
                    (
                        i for i, line in enumerate(result_lines)
                        if i > detailed_cat_start and ('Cis reads' in line or 'Trans reads' in line)
                    ),
                    -1
                )
                if detailed_cat_end != -1 and 'Cis reads' in result_lines[detailed_cat_end]:
                    detailed_cat_end = next(
                        (i for i, line in enumerate(result_lines) if i > detailed_cat_end and 'Trans reads' in line),
                        -1
                    )
                    if detailed_cat_end != -1:
                        detailed_cat_end += 1
                elif detailed_cat_end != -1:
                    detailed_cat_end += 1

                end_idx = detailed_cat_end if detailed_cat_end != -1 else len(result_lines)
                detailed_cat_lines = result_lines[detailed_cat_start:end_idx]
                del result_lines[detailed_cat_start:end_idx]

                qc_status_idx = next(
                    (
                        i for i, line in enumerate(qc_lines)
                        if 'QC PASSED' in line or 'QC passed' in line or 'QC FAILED' in line or 'QC failed' in line
                    ),
                    -1
                )
                insert_idx = qc_status_idx if qc_status_idx != -1 else 0
                qc_lines[insert_idx:insert_idx] = detailed_cat_lines

            chimeric_idx = next((i for i, line in enumerate(result_lines) if 'Chimeric reads percentage:' in line), -1)
            if chimeric_idx != -1:
                chimeric_line = result_lines.pop(chimeric_idx)
                qc_status_idx = next(
                    (
                        i for i, line in enumerate(qc_lines)
                        if 'QC PASSED' in line or 'QC passed' in line or 'QC FAILED' in line or 'QC failed' in line
                    ),
                    -1
                )
                insert_idx = qc_status_idx if qc_status_idx != -1 else 0
                qc_lines[insert_idx:insert_idx] = ['', chimeric_line, '']

        idx_vv = next((i for i, s in enumerate(section_objs) if s['titleNorm'] == 'variant validation'), -1)
        idx_qc = next((i for i, s in enumerate(section_objs) if s['titleNorm'] == 'quality control'), -1)
        if idx_vv != -1:
            if idx_qc == -1:
                section_objs.append({
                    'header': 'Quality control',
                    'title': 'Quality control',
                    'titleNorm': 'quality control',
                    'lines': []
                })
                idx_qc = len(section_objs) - 1
            section_objs[idx_qc]['lines'] = section_objs[idx_vv]['lines'] + section_objs[idx_qc]['lines']
            section_objs.pop(idx_vv)

        rank_map = {
            'result': 0,
            'results': 0,
            'quality control': 1,
            'quality control details': 2,
            'variant calling': 3
        }
        section_objs.sort(key=lambda section: rank_map.get(section['titleNorm'], 100))

        is_phasing_report = any(
            any('Counting' in line and 'WhatsHap' in line and 'HapCUT2' in line for line in section['lines'])
            for section in section_objs
        )
        has_result = any(section['titleNorm'] in ('result', 'results') for section in section_objs)

        if not is_phasing_report and not has_result:
            section_objs.insert(0, {
                'header': 'Result',
                'title': 'Result',
                'titleNorm': 'result',
                'lines': first_content_lines[:]
            })
            section_objs.sort(key=lambda section: rank_map.get(section['titleNorm'], 100))

        idx_result_for_variant_move = next((i for i, s in enumerate(section_objs) if s['titleNorm'] in ('result', 'results')), -1)
        idx_variant_calling = next((i for i, s in enumerate(section_objs) if s['titleNorm'] == 'variant calling'), -1)

        if idx_result_for_variant_move != -1 and idx_variant_calling != -1:
            result_lines = section_objs[idx_result_for_variant_move]['lines']
            variant_lines = section_objs[idx_variant_calling]['lines']
            variant_matching_results_idx = next(
                (i for i, line in enumerate(variant_lines) if 'Variant Matching Results' in line),
                -1
            )
            variant_matching_line = next(
                (
                    line for i, line in enumerate(variant_lines)
                    if 'Variant Matching:' in line
                    and not (
                        variant_matching_results_idx != -1
                        and variant_matching_results_idx < i <= variant_matching_results_idx + 2
                    )
                ),
                ''
            )
            variant1_line = next((line for line in variant_lines if 'Variant 1:' in line), '')

            if variant_matching_line or variant1_line:
                passing_qc_idx = next((i for i, line in enumerate(result_lines) if 'Passing QC reads:' in line), -1)
                insertion = ['']
                if variant_matching_line:
                    insertion.append(variant_matching_line)
                if variant1_line:
                    insertion.append(variant1_line)
                insertion.append('')

                if passing_qc_idx != -1:
                    result_lines[passing_qc_idx + 1:passing_qc_idx + 1] = insertion
                else:
                    result_lines.extend(insertion)

        for section in section_objs:
            if section['titleNorm'] in ('result', 'results'):
                content_text = '\n'.join(section['lines'])
                if 'Counting' in content_text and 'WhatsHap' in content_text and 'HapCUT2' in content_text:
                    section['lines'] = first_content_lines + [''] + section['lines']
                elif not any(line.strip() for line in section['lines']):
                    section['lines'] = first_content_lines[:]
                break

        wanted_titles = {'result', 'results', 'quality control', 'quality control details', 'variant calling'}
        return [section for section in section_objs if section['titleNorm'] in wanted_titles]

    def _line_to_pdf_paragraph(self, line, styles, suffix=''):
        text = escape(line).replace('&lt;-', '<-')
        text = re.sub(r'\b(Cis|Trans)\b', r'<b>\1</b>', text)

        if line.startswith('*') and 'PASSED' in line:
            style = styles[f'PassLine{suffix}']
        elif line.startswith('*'):
            style = styles[f'FailLine{suffix}']
        elif line.endswith('<-'):
            style = styles[f'GreenLine{suffix}']
        else:
            style = styles[f'BodyLine{suffix}']

        return Paragraph(text, style)

    def _generate_sample_report_pdf(self, report_data):
        buffer = BytesIO()
        left_margin = 0.55 * inch
        right_margin = 0.55 * inch
        top_margin = 0.75 * inch
        bottom_margin = 0.70 * inch

        doc = SimpleDocTemplate(
            buffer,
            pagesize=A4,
            rightMargin=right_margin,
            leftMargin=left_margin,
            topMargin=top_margin,
            bottomMargin=bottom_margin
        )
        base_styles = getSampleStyleSheet()
        styles = {
            'SampleIDHeading': ParagraphStyle(
                'SampleIDHeading',
                parent=base_styles['Normal'],
                fontSize=13,
                leading=17,
                alignment=TA_RIGHT,
                spaceAfter=2
            ),
            'AnalysisDate': ParagraphStyle(
                'AnalysisDate',
                parent=base_styles['Normal'],
                fontSize=9,
                leading=12,
                alignment=TA_RIGHT,
                textColor=colors.HexColor('#666666'),
                spaceAfter=10
            ),
            'SectionTitle': ParagraphStyle(
                'SectionTitle',
                parent=base_styles['Heading2'],
                textColor=colors.HexColor('#0b4b63'),
                fontSize=14,
                leading=18,
                spaceBefore=14,
                spaceAfter=8
            ),
            'BodyLine': ParagraphStyle(
                'BodyLine',
                parent=base_styles['BodyText'],
                fontSize=9,
                leading=12,
                spaceAfter=2
            ),
            'PassLine': ParagraphStyle(
                'PassLine',
                parent=base_styles['BodyText'],
                textColor=colors.HexColor('#017a01'),
                fontSize=9,
                leading=12,
                spaceAfter=2
            ),
            'FailLine': ParagraphStyle(
                'FailLine',
                parent=base_styles['BodyText'],
                textColor=colors.red,
                fontSize=9,
                leading=12,
                spaceAfter=2
            ),
            'GreenLine': ParagraphStyle(
                'GreenLine',
                parent=base_styles['BodyText'],
                textColor=colors.HexColor('#1d8214'),
                fontSize=9,
                leading=12,
                spaceAfter=2
            ),
            'BodyLineLg': ParagraphStyle(
                'BodyLineLg',
                parent=base_styles['BodyText'],
                fontSize=10,
                leading=13,
                spaceAfter=2
            ),
            'PassLineLg': ParagraphStyle(
                'PassLineLg',
                parent=base_styles['BodyText'],
                textColor=colors.HexColor('#017a01'),
                fontSize=10,
                leading=13,
                spaceAfter=2
            ),
            'FailLineLg': ParagraphStyle(
                'FailLineLg',
                parent=base_styles['BodyText'],
                textColor=colors.red,
                fontSize=10,
                leading=13,
                spaceAfter=2
            ),
            'GreenLineLg': ParagraphStyle(
                'GreenLineLg',
                parent=base_styles['BodyText'],
                textColor=colors.HexColor('#1d8214'),
                fontSize=10,
                leading=13,
                spaceAfter=2
            )
        }

        ampmap_version = report_data['ampmap_version']
        sample_name = report_data['displayed_sample_name']
        analysis_datetime = report_data.get('analysis_datetime', '')
        page_width, page_height = A4
        gray = colors.HexColor('#888888')
        line_gray = colors.HexColor('#cccccc')

        class ReportCanvas(pdfcanvas.Canvas):
            def __init__(self, *args, **kwargs):
                pdfcanvas.Canvas.__init__(self, *args, **kwargs)
                self._saved_page_states = []

            def showPage(self):
                self._saved_page_states.append(dict(self.__dict__))
                self._startPage()

            def save(self):
                num_pages = len(self._saved_page_states)
                for state in self._saved_page_states:
                    self.__dict__.update(state)
                    self._draw_header_footer(num_pages)
                    pdfcanvas.Canvas.showPage(self)
                pdfcanvas.Canvas.save(self)

            def _draw_header_footer(self, page_count):
                self.saveState()

                # Header
                header_y = page_height - 0.42 * inch
                self.setFont('Helvetica', 8)
                self.setFillColor(gray)
                self.drawString(left_margin, header_y, f"{ampmap_version} Report")
                self.drawRightString(page_width - right_margin, header_y, f"Sample ID: {sample_name}")
                self.setStrokeColor(line_gray)
                self.setLineWidth(0.5)
                self.line(left_margin, header_y - 5, page_width - right_margin, header_y - 5)

                # Footer
                footer_y = 0.38 * inch
                self.setFont('Helvetica', 8)
                self.setFillColor(gray)
                self.drawString(left_margin, footer_y, "NSW Health Pathology Randwick Genetics")
                self.drawRightString(
                    page_width - right_margin, footer_y,
                    f"Page {self._pageNumber} of {page_count}"
                )
                self.setStrokeColor(line_gray)
                self.line(left_margin, footer_y + 12, page_width - right_margin, footer_y + 12)

                self.restoreState()

        story = [
            Paragraph(escape(sample_name), styles['SampleIDHeading']),
            Paragraph(f"Analysis date: {escape(analysis_datetime)}", styles['AnalysisDate']),
            Spacer(1, 0.10 * inch)
        ]

        for section in self._build_sample_report_section_objects(report_data['report_content']):
            title_norm = section['titleNorm']
            is_large_section = title_norm in ('result', 'results', 'quality control')
            suffix = 'Lg' if is_large_section else ''

            story.append(Paragraph(escape(section['title']), styles['SectionTitle']))
            lines = section['lines'] or ['No data available.']
            for line in lines:
                if title_norm in ('result', 'results') and line.strip().startswith('Variant'):
                    continue
                if not line.strip():
                    story.append(Spacer(1, 0.08 * inch))
                else:
                    story.append(self._line_to_pdf_paragraph(line, styles, suffix))

        doc.build(story, canvasmaker=ReportCanvas)
        buffer.seek(0)
        return buffer

    def init_routes(self):
        """Initialize Flask routes"""
        
        @self.app.route('/')
        def index():
            return render_template('index.html')

        @self.app.route('/api/runs')
        def get_runs():
            try:
                runs = [d for d in os.listdir(str(LOCAL_PATH)) 
                        if os.path.isdir(os.path.join(str(LOCAL_PATH), d))]
                return jsonify(sorted(runs))
            except Exception as e:
                return jsonify({'error': str(e)}), 500

        @self.app.route('/api/runs/<run_id>/html')
        def get_run_html(run_id):
            try:
                run_dir = os.path.join(str(LOCAL_PATH), run_id)
                html_files = glob.glob(os.path.join(run_dir, '*.html'))
                if not html_files:
                    return jsonify({'error': 'HTML file not found'}), 404
                
                with open(html_files[0], 'r') as f:
                    content = f.read()
                return content
            except Exception as e:
                return jsonify({'error': str(e)}), 500

        @self.app.route('/api/runs/<run_id>/samples')
        def get_samples(run_id):
            try:
                samples = set()
                result_dir = os.path.join(str(LOCAL_PATH), run_id, 'result')
                
                for i in range(1, 93):
                    barcode_dir = os.path.join(result_dir, f'barcode{i:02d}')
                    if not os.path.exists(barcode_dir):
                        continue
                        
                    report_files = glob.glob(os.path.join(barcode_dir, '*_report.txt'))
                    for report_file in report_files:
                        episode = os.path.basename(report_file).replace('_report.txt', '')
                        if episode:
                            samples.add(f"{episode}_b{i:02d}")
                
                return jsonify(sorted(list(samples)))
            except Exception as e:
                return jsonify({'error': str(e)}), 500

        @self.app.route('/api/runs/<run_id>/samples/<sample_id>/report')
        def get_sample_report(run_id, sample_id):
            try:
                report_data, status_code = self._get_sample_report_data(run_id, sample_id)
                if status_code != 200:
                    return jsonify(report_data), status_code

                return report_data['report_content']
            except Exception as e:
                return jsonify({'error': str(e)}), 500

        @self.app.route('/api/runs/<run_id>/samples/<sample_id>/download/<file_type>')
        def download_file(run_id, sample_id, file_type):
            try:
                result_dir = os.path.join(str(LOCAL_PATH), run_id, 'result')
                
                match = re.match(r"(.+)_b(\d{2})$", sample_id)
                if not match:
                    return jsonify({'error': 'Invalid sample ID format'}), 400
                
                original_sample_id, barcode = match.groups()
                
                file_names = {
                    'xml': f'{original_sample_id}.xml',
                    'HapCUT2': f'HapCUT2.log',
                    'WhatsHap': f'whatshap.log',
                    'Pipeline': f'pipeline.log',
                    'MainPipeline': f'{run_id}.log',
                    'report': f'{original_sample_id}_AmpMap_report.pdf'
                }
                
                if file_type not in file_names:
                    return jsonify({'error': 'Invalid file type'}), 400
                    
                file_name = file_names[file_type]

                if file_type == 'report':
                    report_data, status_code = self._get_sample_report_data(run_id, sample_id)
                    if status_code != 200:
                        return jsonify(report_data), status_code

                    pdf_buffer = self._generate_sample_report_pdf(report_data)
                    safe_file_name = re.sub(r'[^\w.\-]+', '_', file_name)
                    return send_file(
                        pdf_buffer,
                        as_attachment=True,
                        download_name=safe_file_name,
                        mimetype='application/pdf'
                    )
                
                if file_type == 'MainPipeline':
                    # Main pipeline log is in result directory
                    file_path = os.path.join(result_dir, file_name)
                else:
                    # Other files are in barcode directory
                    barcode_dir = os.path.join(result_dir, f'barcode{barcode}')
                    file_path = os.path.join(barcode_dir, file_name)
                
                if not os.path.exists(file_path):
                    return jsonify({'error': 'File not found'}), 404
                    
                return send_file(file_path, as_attachment=True)
            except Exception as e:
                return jsonify({'error': str(e)}), 500



    def start(self):
        """Start both the monitoring thread and Flask app"""
        try:
            self.logger.info("Starting GUI AmpMap components")
            
            # Start remote monitoring in a separate thread
            monitoring_thread = Thread(target=self.remote_monitor.run)
            monitoring_thread.daemon = True
            monitoring_thread.start()
            self.logger.info("Remote monitoring thread started")
            
            # Start Flask app
            self.logger.info(f"Starting web interface on {FLASK_HOST}:{FLASK_PORT}")
            self.app.run(host=FLASK_HOST, port=FLASK_PORT)
            
        except Exception as e:
            self.logger.error(f"Error starting GUI AmpMap: {str(e)}")
            raise

if __name__ == '__main__':
    # Set up logging configuration
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    
    # Ensure log directory exists
    log_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Create logs directory
    logs_dir = os.path.join(log_dir, 'logs')
    os.makedirs(logs_dir, exist_ok=True)
    
    # Set up common file handler for all logs
    file_handler = logging.FileHandler(os.path.join(logs_dir, 'gui_ampmap.log'))
    file_handler.setFormatter(logging.Formatter(log_format))
    
    # Configure root logger
    logging.basicConfig(
        level=logging.INFO,
        format=log_format,
        handlers=[
            logging.StreamHandler(),  # Console output
            file_handler  # File output
        ]
    )
    
    # Configure component-specific loggers
    remote_monitor_logger = logging.getLogger('RemoteMonitor')
    remote_monitor_logger.setLevel(logging.INFO)
    remote_monitor_logger.propagate = False  # Prevent propagation to root logger
    remote_monitor_handler = logging.FileHandler(os.path.join(logs_dir, 'remotemonitor.log'))
    remote_monitor_handler.setFormatter(logging.Formatter(log_format))
    remote_monitor_logger.addHandler(remote_monitor_handler)
    
    # GUI_AmpMap logger - separate from root logger
    gui_logger = logging.getLogger('GUI_AmpMap')
    gui_logger.setLevel(logging.INFO)
    gui_logger.propagate = False  # Prevent propagation to root logger
    gui_handler = logging.FileHandler(os.path.join(logs_dir, 'gui_ampmap.log'))
    gui_handler.setFormatter(logging.Formatter(log_format))
    gui_logger.addHandler(gui_handler)
    gui_logger.addHandler(logging.StreamHandler())  # Console output for GUI
    
    root_logger = logging.getLogger()
    root_logger.info("=== GUI AmpMap Starting ===")
    root_logger.info(f"Log directory: {logs_dir}")
    
    try:
        # Start the application
        app = GUI_AmpMap()
        app.start()
    except Exception as e:
        root_logger.error(f"Fatal error: {str(e)}", exc_info=True)
        raise
