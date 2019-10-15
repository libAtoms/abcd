import logging
from abcd.server.app import create_app

logging.basicConfig(level=logging.DEBUG)
app = create_app()

if __name__ == '__main__':
    # app.logger.setLevel(logging.DEBUG)
    # handler = logging.StreamHandler()
    # app.logger.addHandler(handler)

    app.run(host='0.0.0.0', debug=True)
