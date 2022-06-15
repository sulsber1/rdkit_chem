FROM python:3.8

WORKDIR /rdkit_api
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY . .