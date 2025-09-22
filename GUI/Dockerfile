FROM python:3.9-slim

WORKDIR /app

# Copy requirements and install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application files
COPY . .

# Create directories for data and logs
RUN mkdir -p /app/data /app/logs

# Expose port
EXPOSE 5001

# Run the application
CMD ["python", "gui_ampmap.py"]