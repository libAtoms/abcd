from abcd_server.app import create_app

app = create_app()

if __name__ == '__main__':
    print('hello')
    app.run()
