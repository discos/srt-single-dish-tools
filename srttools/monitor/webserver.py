import asyncio
import base64
import json
import threading
import warnings
import re
import os
from string import Template

try:
    from importlib.resources import files
except ImportError:
    from importlib_resources import files

try:
    import tornado.web
    import tornado.websocket
    from tornado.web import Application, RequestHandler
    from tornado.websocket import WebSocketHandler
except ImportError:
    raise ImportError(
        "To use SDTmonitor, you need to install tornado: \n" "\n   > pip install tornado"
    )

from srttools.monitor.common import log


def create_index_file(**kwargs):
    config = {
        "port": 8080,
        "reconnectTime": 10000,
        "maxScrollSpeed": 80,
    }
    config.update(kwargs)
    with open("index.html", "w") as webpage:
        with open(files("srttools.monitor").joinpath("resources", "index.html"), "r") as template:
            print(Template(template.read()).safe_substitute(config), file=webpage, end="")


class WSHandler(WebSocketHandler):
    # This allows clients that did not send any request to the HTTPHandler previously
    # i.e.: a client that opens the index.html page instead of accessing it via network
    check_origin = lambda _, __: True
    # We don't care about listening to messages
    on_message = lambda _, __: None

    def initialize(self, connected_clients, images):
        self.connected_clients = connected_clients
        self.images = images

    def open(self):
        log.info(f"Got connection from {self.request.remote_ip}")
        self.connected_clients.add(self)
        # Send all the images to new clients
        keys = self.images.keys()
        for index in keys:
            self.send_image(index)

    def on_close(self):
        self._close()

    def send_image(self, index):
        message = {"index": index, "image": self.images.get(index, "")}
        try:
            self.write_message(json.dumps(message))
        except (
            tornado.websocket.WebSocketClosedError,
            tornado.iostream.StreamClosedError,
        ):  # pragma: no cover
            # Tried to cover this exception by dropping the connection in multiple ways after it was opened but it was never covered, ignoring it
            self._close()

    def _close(self):
        if self in self.connected_clients:
            self.connected_clients.remove(self)
            log.info(f"Client {self.request.remote_ip} disconnected")


class HTTPHandler(RequestHandler):
    def get(self):
        # Answer the HTTP request with the index.html page
        self.write(open("index.html").read())


class FaviconHandler(RequestHandler):
    def initialize(self):
        # Load the favicon in memory
        favicon_data = None
        with files("srttools.monitor").joinpath("resources", "favicon.ico").open("rb") as f:
            favicon_data = f.read()
        self.favicon_data = favicon_data

    def get(self):
        self.set_header("Content-Type", "image/x-icon")
        self.write(self.favicon_data)


class WebServer:
    def __init__(self, extension, localhost=False, port=8080):
        self.extension = extension
        self.port = port
        self.address = "127.0.0.1" if localhost else "0.0.0.0"

        # Load the current images
        self.images = {}
        pattern = re.compile(rf"^latest_\d{{3}}\.{extension}$")
        files = [f for f in os.listdir() if pattern.match(f)]
        files.sort()
        for filename in files:
            self._load_image(filename)

        self.connected_clients = set()

        self.t = None
        self.started = False
        application = Application(
            [
                (
                    r"/images",
                    WSHandler,
                    dict(
                        connected_clients=self.connected_clients,
                        images=self.images,
                    ),
                ),
                (r"/", HTTPHandler),
                (r"/index.html", HTTPHandler),
                (r"/favicon.ico", FaviconHandler),
            ]
        )

        # Disable default log function, we use custom ones
        def log_function(_):
            pass

        application.log_request = log_function
        self.web_server = tornado.httpserver.HTTPServer(application)
        asyncio.set_event_loop(asyncio.new_event_loop())
        try:
            self.web_server.listen(self.port, address=self.address)
        except OSError:
            raise OSError(f"Port {self.port} is already being used, choose a different one!")

    def start(self):
        create_index_file(port=self.port)

        self.ioloop = tornado.ioloop.IOLoop.current()
        if not self.started and not self.ioloop.asyncio_loop.is_running():
            self.t = threading.Thread(target=self.ioloop.start)
            self.t.start()
        self.started = True

    def stop(self):
        if self.started:
            self.ioloop.add_callback(self.ioloop.stop)
            if self.t:
                self.t.join()
        self.started = False

    def _load_image(self, image_file):
        index = int(image_file.split("_")[1].split(".")[0])
        try:
            image_string = base64.b64encode(open(image_file, "rb").read())
            image_string = image_string.decode("utf-8")
            self.images[index] = image_string
        except OSError:
            self.images.pop(index, None)
            image_string = ""
        return index, image_string

    def update(self, image_file):
        # Update the image in memory before sending it
        index, image_string = self._load_image(image_file)
        clients = self.connected_clients
        for client in clients:
            self.ioloop.add_callback(client.send_image, index)
