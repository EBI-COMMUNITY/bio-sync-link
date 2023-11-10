FROM python:3.11
LABEL authors="mernst"

VOLUME /app

WORKDIR /app

COPY . /app

RUN pip install -r requirements.txt

EXPOSE 80

CMD ["python", "run.py"]
