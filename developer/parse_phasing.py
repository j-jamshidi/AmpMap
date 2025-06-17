#!/usr/bin/env python3
import os
import re
from pathlib import Path

def extract_value(line, field_type):
    """Extract values based on field type."""
    if field_type == "variant":
        # Extract text after "-"
        if "-" in line:
            return line.split("-", 1)[1].strip()
        return ""
    
    elif field_type == "phase_determination":
        # Extract the last word
        words = line.strip().split()
        if words:
            return words[-1]
        return ""
    
    elif field_type == "number_with_percentage":
        # Extract both number and percentage
        number_match = re.search(r'(\d+,?\d*)\s*\(([\d.]+)%\)', line)
        if number_match:
            return (number_match.group(1).replace(',', ''), number_match.group(2))
        return ("", "")
    
    else:  # simple number
        number_match = re.search(r'(\d+,?\d*)', line)
        if number_match:
            return number_match.group(1).replace(',', '')
        return ""

def parse_report_file(filepath):
    """Parse a single report file and extract specified fields."""
    data = {}
    field_mapping = {
        "Variant 1:": ("Variant 1", "variant"),
        "Variant 2:": ("Variant 2", "variant"),
        "Number of variants called from the amplicon:": ("Number of variants", "number"),
        "Total reads:": ("Total reads", "number"),
        "Total high quality spanning reads:": ("Total high quality spanning reads", "number"),
        "Reads with ref allele for both variants (Cis):": ("Reads with ref allele for both variants", "number_with_percentage"),
        "Reads with alt allele for both variants (Cis):": ("Reads with alt allele for both variants", "number_with_percentage"),
        "Reads with ref for first, alt for second (Trans):": ("Reads with ref first alt second", "number_with_percentage"),
        "Reads with alt for first, ref for second (Trans):": ("Reads with alt first ref second", "number_with_percentage"),
        "Cis reads (both ref or both alt):": ("Cis reads", "number_with_percentage"),
        "Trans reads (one ref, one alt):": ("Trans reads", "number_with_percentage"),
        "Chimeric reads percentage:": ("Chimeric reads percentage", "number"),
        "Counting": ("Counting", "phase_determination"),
        "WhatsHap": ("WhatsHap", "phase_determination"),
        "HapCUT2": ("HapCUT2", "phase_determination")
    }
    
    with open(filepath, 'r') as f:
        for line in f:
            for search_text, (header, field_type) in field_mapping.items():
                if line.strip().startswith(search_text):
                    value = extract_value(line, field_type)
                    if field_type == "number_with_percentage":
                        data[header] = value[0]  # number
                        data[f"{header} %"] = value[1]  # percentage
                    else:
                        data[header] = value
                    break
    
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
    
    # Define the order of columns
    columns = [
        "Variant 1",
        "Variant 2",
        "Number of variants",
        "Total reads",
        "Total high quality spanning reads",
        "Reads with ref allele for both variants",
        "Reads with ref allele for both variants %",
        "Reads with alt allele for both variants",
        "Reads with alt allele for both variants %",
        "Reads with ref first alt second",
        "Reads with ref first alt second %",
        "Reads with alt first ref second",
        "Reads with alt first ref second %",
        "Cis reads",
        "Cis reads %",
        "Trans reads",
        "Trans reads %",
        "Chimeric reads percentage",
        "Counting",
        "WhatsHap",
        "HapCUT2"
    ]
    
    # Parse all files
    all_data = []
    for filepath in report_files:
        data = parse_report_file(filepath)
        barcode = filepath.parent.name
        episode = filepath.stem.replace('_report', '')
        all_data.append((barcode, episode, data))
    
    # Write output file
    with open('parsed_phasing_reports.tsv', 'w') as f:
        # Write header line
        header_line = ['barcode', 'episode'] + columns
        f.write('\t'.join(header_line) + '\n')
        
        # Write data lines
        for barcode, episode, data in all_data:
            line_values = [barcode, episode]
            for column in columns:
                line_values.append(str(data.get(column, '')))
            f.write('\t'.join(line_values) + '\n')

if __name__ == '__main__':
    main()