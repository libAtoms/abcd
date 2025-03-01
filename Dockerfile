FROM python:3
#FROM python:3.7.2-alpine

COPY abcd/server/requirements.txt /
RUN pip install --no-cache-dir -r /requirements.txt

COPY . /usr/src/app
WORKDIR /usr/src/app

RUN python setup.py develop

#EXPOSE 5000
#CMD ["/bin/bash"]
CMD [ "python", "abcd/server/__init__.py" ]
#CMD ["gunicorn", "-w 4", "-b 0.0.0.0:8000", "server:app"]
