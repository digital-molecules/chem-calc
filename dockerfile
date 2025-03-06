# Use the official Streamlit base image or Python base image
FROM python:3.12-slim

# Install necessary system dependencies for graphical libraries
RUN apt-get update && apt-get install -y \
    libxrender1 \
    libx11-dev \
    libxext-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
COPY requirements.txt .
RUN pip install -r requirements.txt

# Copy the app code into the container
COPY . /app
WORKDIR /app

# Expose the app's port
EXPOSE 8501

# Run the Streamlit app
CMD ["streamlit", "run", "main.py"]
