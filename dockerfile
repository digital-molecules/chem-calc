FROM python:3.12-slim

RUN apt-get update && apt-get install -y \
    libxrender1 \
    libx11-dev \
    libxext-dev \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN pip install -r requirements.txt

COPY . /app
WORKDIR /app

EXPOSE 8501

CMD ["streamlit", "run", "main.py"]
