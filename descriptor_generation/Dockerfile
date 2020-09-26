FROM continuumio/miniconda3:4.8.2

RUN /opt/conda/bin/conda install -c rdkit -c mordred-descriptor mordred -y --quiet
RUN python -m pip install Flask

COPY ./app/ /app

CMD python /app/app.py