import threading
import warnings
import asyncio
import json
import base64

try:
    import tornado.web
    import tornado.websocket
    from tornado.web import Application, RequestHandler
    from tornado.websocket import WebSocketHandler
except ImportError:
    warnings.warn("To use SDTmonitor, you need to install tornado: \n" "\n   > pip install tornado")
    RequestHandler = WebSocketHandler = object

from srttools.monitor.common import MAX_FEEDS, log


def create_index_file(port, max_images=MAX_FEEDS * 2):
    html_string = (
        """
<!DOCTYPE html>
<html>
    <head>
        <meta charset="UTF-8">
        <title>SRT Quicklook</title>
    </head>
    <body>
        <script type="text/javascript">
            window.onload = function()
            {
                function init_images(n)
                {
                    if(n % 2 == 0)
                    {
                        n++;
                    }

                    for(i = 0; i <= n; i++)
                    {
                        var div = document.getElementById("div_" + i.toString());

                        if(div == null)
                        {
                            var image = new Image();
                            image.id = "image_" + i.toString();
                            image.style.width = "100%";
                            image.src = "data:image/gif;base64,R0lGODlhAQABAAD/ACwAAAAAAQABAAACADs%3D";

                            div = document.createElement("DIV");
                            div.setAttribute("id", "div_" + i.toString());
                            div.setAttribute("style", "width:50%; float:left;");

                            div.appendChild(image);
                            document.body.appendChild(div);
                        }
                    }
                }

                function set_visibility()
                {
                    var white_image = "data:image/gif;base64,R0lGODlhAQABAAD/ACwAAAAAAQABAAACADs%3D";

                    for(i = 0; i < document.getElementsByTagName("IMG").length; i+=2)
                    {
                        var left = i;
                        var right = i+1;
                        var left_div = document.getElementById("div_" + left.toString());
                        var right_div = document.getElementById("div_" + right.toString());
                        var left_image = document.getElementById("image_" + left.toString());
                        var right_image = document.getElementById("image_" + right.toString());

                        if(left_image.src == white_image && right_image.src == white_image)
                        {
                            left_image.style.display = "none";
                            right_image.style.display = "none";
                        }
                        else
                        {
                            left_image.style.display = "block";
                            right_image.style.display = "block";
                        }
                    }
                }

                function connect()
                {
                    var destination = document.location.href;
                    if(destination.startsWith("file"))
                    {
                        destination = "localhost";
                    }
                    else
                    {
                        destination = document.location.href.split(":")[1]
                    }

                    var ws = new WebSocket("ws:" + destination + ":"""
        + str(port)
        + """/images");

                    ws.onopen = function()
                    {
                        console.log('Connected')
                    }

                    ws.onmessage = function(message)
                    {
                        var msg = JSON.parse(message.data)
                        init_images(msg.index);

                        var image = document.getElementById("image_" + msg.index.toString());

                        if(msg.image == "")
                        {
                            image.src = "data:image/gif;base64,R0lGODlhAQABAAD/ACwAAAAAAQABAAACADs%3D";
                        }
                        else
                        {
                            image.src = "data:image/png;base64," + msg.image;
                        }

                        set_visibility();
                    };

                    ws.onclose = function(e)
                    {
                        console.log('Socket is closed. Reconnect will be attempted in 10 seconds.');
                        setTimeout(function()
                        {
                            connect();
                        }, 10000);
                    };

                    ws.onerror = function(err)
                    {
                        console.error('Socket encountered error. Closing socket');
                        ws.close();
                    };
                }

                connect();
            }
        </script>
    </body>
</html>"""
    )
    with open("index.html", "w") as fobj:
        print(html_string, file=fobj)


class WSHandler(WebSocketHandler):
    def initialize(self, connected_clients, images):
        self.connected_clients = connected_clients
        self.images = images

    def check_origin(self, origin):
        # This allows clients that did not send any request to the HTTPHandler previously
        # i.e.: a client that opens the index.html page instead of accessing it via network
        return True

    def open(self):
        log.info("Got connection from {}".format(self.request.remote_ip))
        self.connected_clients.add(self)
        # Send all the images to new clients
        keys = self.images.keys()
        for index in keys:
            self.send_image(index)

    def on_close(self):
        self._close()

    def on_message(self, message):
        pass

    def send_image(self, index):
        message = {"index": index, "image": self.images[index]}
        try:
            self.write_message(json.dumps(message))
        except tornado.websocket.WebSocketClosedError:
            self._close()

    def _close(self):
        if self in self.connected_clients:
            self.connected_clients.remove(self)
            log.info("Client {} disconnected".format(self.request.remote_ip))


class HTTPHandler(RequestHandler):
    def get(self):
        # Answer the HTTP request with the index.html page
        self.write(open("index.html", "r").read())


class WebServer(object):
    def __init__(self, extension, port=8080):
        self.extension = extension
        self.port = port

        # Load the current images
        self.images = {}
        for index in range(MAX_FEEDS * 2):
            self._load_image("latest_{}.{}".format(index, extension))

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
            ]
        )

        # Disable default log function, we use custom ones
        def log_function(_):
            pass

        application.log_request = log_function
        self.web_server = tornado.httpserver.HTTPServer(application)
        try:
            self.web_server.listen(self.port)
        except OSError:
            raise OSError(
                "Port {} is already being used, choose a different one!".format(self.port)
            )

    def start(self):
        self._asyncioloop = None
        try:
            asyncio.get_event_loop()
        except RuntimeError:
            self._asyncioloop = asyncio.new_event_loop()
            asyncio.set_event_loop(self._asyncioloop)
        self.ioloop = tornado.ioloop.IOLoop.current()

        create_index_file(self.port)

        self.t = threading.Thread(target=self.ioloop.start)
        self.t.start()
        self.started = True

    def stop(self):
        if self.started:
            self.ioloop.add_callback(self.ioloop.stop)
            if self._asyncioloop:
                self._asyncioloop.stop()
        if self.t:
            self.t.join()
        self.started = False

    def _load_image(self, image_file):
        index = int(image_file.split("_")[1].split(".")[0])
        try:
            image_string = base64.b64encode(open(image_file, "rb").read())
            image_string = image_string.decode("utf-8")
        except IOError:
            image_string = ""
        self.images[index] = image_string
        return index, image_string

    def update(self, image_file):
        # Update the image in memory before sending it
        index, image_string = self._load_image(image_file)
        clients = self.connected_clients
        for client in clients:
            self.ioloop.add_callback(client.send_image, index)
