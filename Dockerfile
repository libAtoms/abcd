FROM python:3

COPY abcd_server/requirements.txt /
RUN pip install --no-cache-dir -r /requirements.txt

COPY . /usr/src/app
WORKDIR /usr/src/app

RUN python setup.py develop

#EXPOSE 5000
#CMD ["/bin/bash"]
#CMD [ "python", "./your-daemon-or-script.py" ]
CMD ["gunicorn", "-w 4", "-b 0.0.0.0:8000", "abcd_server:app"]


