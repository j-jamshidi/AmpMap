#!/usr/bin/env python3
import os
import re
from pathlib import Path

def clean_header(header):
    """Clean header by removing special characters and spaces."""
    return header.strip().replace('≥', 'ge')

def parse_report_file(filepath):
    """Parse a single report file and return a dictionary of values."""
    data = {}
    with open(filepath, 'r') as f:
        current_section = None
        for line in f:
            line = line.strip()
            if line and ':' in line:
                # Get header and value
                header, value = [x.strip() for x in line.split(':', 1)]
                
                # Special handling for Amplicon region
                if header == "Amplicon region":
                    data[header] = value
                    continue
                    
                # Special handling for Amplicon length
                if header == "Amplicon length":
                    bp_match = re.search(r'(\d+,?\d*)', value)
                    if bp_match:
                        data[header] = bp_match.group(1).replace(',', '')
                    continue
                
                # Handle parentheses percentages
                parentheses_match = re.search(r'(.*?)\s*\(([\d.]+)%\)', value)
                if parentheses_match:
                    main_value, percentage = parentheses_match.groups()
                    data[header] = main_value.replace(',', '')
                    data[f"{header} %"] = percentage
                    continue
                
                # Remove any unit indicators (X, bp) and commas
                value = value.replace(',', '').replace('X', '').replace('bp', '').strip()
                data[header] = value
                
    return data

def main():
    # Get all report files
    report_files = []
    for path in Path('.').rglob('*_report.txt'):
        if path.parent.name.startswith('barcode'):
            report_files.append(path)
    
    if not report_files:
        print("No report files found!")
        return

    # Parse first file to get all possible headers
    headers = set()
    first_data = parse_report_file(report_files[0])
    headers.update(first_data.keys())
    
    # Parse all files to ensure we have all possible headers
    all_data = []
    for filepath in report_files:
        data = parse_report_file(filepath)
        headers.update(data.keys())
        barcode = filepath.parent.name
        episode = filepath.stem.replace('_report', '')
        all_data.append((barcode, episode, data))
    
    # Sort headers for consistent output
    headers = sorted(headers)
    
    # Write output file
    with open('parsed_reports.tsv', 'w') as f:
        # Write header line
        header_line = ['barcode', 'episode'] + [clean_header(h) for h in headers]
        f.write('\t'.join(header_line) + '\n')
        
        # Write data lines
        for barcode, episode, data in all_data:
            line_values = [barcode, episode]
            for header in headers:
                line_values.append(str(data.get(header, '')))
            f.write('\t'.join(line_values) + '\n')

if __name__ == '__main__':
    main()