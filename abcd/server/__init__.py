from abcd.server.app import create_app

app = create_app()

if __name__ == '__main__':
    print('hello')
    app.run(host='0.0.0.0')
