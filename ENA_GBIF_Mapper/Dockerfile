FROM python:3.11
LABEL authors="mernst"

VOLUME /app

WORKDIR /app

COPY . /app

RUN pip install -r requirements.txt

# EXPOSE 80

RUN adduser --disabled-password --gecos '' --system --uid 1001 python && chown -R python /app/run.py

USER 1001

CMD ["python", "run.py"]

