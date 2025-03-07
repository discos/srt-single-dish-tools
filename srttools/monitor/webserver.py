import asyncio
import base64
import json
import threading
import warnings

try:
    import tornado.web
    import tornado.websocket
    from tornado.web import Application, RequestHandler
    from tornado.websocket import WebSocketHandler
except ImportError:
    warnings.warn("To use SDTmonitor, you need to install tornado: \n" "\n   > pip install tornado")
    RequestHandler = WebSocketHandler = object

from srttools.monitor.common import MAX_FEEDS, log


def create_index_file(port):
    html_string = (
        """<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>SRT Quicklook</title>
    <style>
        body {
            margin: 0;
            overflow: hidden;
            background: white;
        }
        #main-view {
            width: 100%;
            height: 90vh;
            display: flex;
            justify-content: center;
            align-items: center;
        }
        .main-image {
            width: auto;
            height: 100%;
            max-width: 50%;
            object-fit: contain;
        }
        #thumbnail-bar {
            width: 100%;
            height: 10vh;
            overflow-x: auto;
            overflow-y: hidden;
            white-space: nowrap;
            background: white;
            display: flex;
            align-items: center;
            justify-content: center;
        }
        .thumbnail-pair {
            display: inline-block;
            height: 10vh;
            line-height: 10vh;
            cursor: pointer;
            object-fit: cover;
            text-align: center;
        }
        .thumbnail {
            z-index: 0;
            position: relative;
            max-height: 10vh;
            object-fit: cover;
        }
        .feed_id {
            z-index: 1;
            position: relative;
            height: 0vh;
            font-family: 'DejaVu Sans', sans-serif;
        }
    </style>
</head>
<body>
    <div id="main-view">
        <img id="main-left" class="main-image" src="" />
        <img id="main-right" class="main-image" src="" />
    </div>
    <div id="thumbnail-bar"></div>

    <script>
        const whiteImage = "data:image/gif;base64,R0lGODlhAQABAAD/ACwAAAAAAQABAAACADs%3D";
        let selectedPairIndex = 0;
        let images = [];
        for (let index = 0; index < """
        + str(MAX_FEEDS)
        + """ * 2; index++) {
            images[index] = whiteImage;
        }
        document.getElementById("main-left").src = whiteImage;
        document.getElementById("main-right").src = whiteImage;

        function changePairsVisibility() {
            let firstAvailable = """
        + str(MAX_FEEDS)
        + """;

            for (let pairIndex = 0; pairIndex < """
        + str(MAX_FEEDS)
        + """; pairIndex++) {
                let pair = document.getElementById("thumb-pair-" + pairIndex);
                if (pair) {
                    let leftIndex = pairIndex * 2;
                    let rightIndex = leftIndex + 1;
                    if (images[leftIndex] === whiteImage && images[rightIndex] === whiteImage) {
                        document.getElementById("thumb-pair-" + pairIndex).hidden = true;
                        document.getElementById("feed-" + pairIndex).innerText = "";
                    } else {
                        document.getElementById("thumb-pair-" + pairIndex).hidden = false;
                        document.getElementById("feed-" + pairIndex).innerText = "FEED " + pairIndex;
                        firstAvailable = pairIndex;
                    }
                }
            }

            if (!document.getElementById("thumb-pair-" + selectedPairIndex).hidden) {
                return selectedPairIndex;
            } else {
                return firstAvailable;
            }
        }

        function addPair(pairIndex) {
            if (!document.getElementById("thumb-pair-" + pairIndex)) {
                let div = document.createElement("div");
                div.classList.add("thumbnail-pair");
                div.id = "thumb-pair-" + pairIndex;
                div.setAttribute("data-pair", pairIndex);
                let leftIndex = pairIndex * 2;
                let rightIndex = leftIndex + 1;
                div.setAttribute("data-left", leftIndex);
                div.setAttribute("data-right", rightIndex);
                div.onclick = function() {
                    selectedPairIndex = parseInt(this.getAttribute("data-pair"));
                    updateMainView();
                };

                let img1 = document.createElement("img");
                img1.classList.add("thumbnail");
                img1.id = "thumb-" + leftIndex;
                img1.src = images[leftIndex];

                let img2 = document.createElement("img");
                img2.classList.add("thumbnail");
                img2.id = "thumb-" + rightIndex;
                img2.src = images[rightIndex];

                let feed_id = document.createElement("div");
                feed_id.classList.add("feed_id");
                feed_id.id = "feed-" + pairIndex;

                div.appendChild(feed_id);
                div.appendChild(img1);
                div.appendChild(img2);
                document.getElementById("thumbnail-bar").appendChild(div);
            }
        }

        function updateMainView() {
            let mainPair = changePairsVisibility();

            if (mainPair !== """
        + str(MAX_FEEDS)
        + """) {
                let leftIndex = mainPair * 2;
                let rightIndex = leftIndex + 1;
                document.getElementById("main-left").src = images[leftIndex];
                document.getElementById("main-right").src = images[rightIndex];
            } else {
                document.getElementById("main-left").src = whiteImage;
                document.getElementById("main-right").src = whiteImage;
            }
        }

        function addImage(index, base64) {
            images[index] = base64 ? "data:image/png;base64," + base64 : whiteImage;

            let pairIndex = Math.floor(index / 2);
            let leftIndex = pairIndex * 2;
            let rightIndex = leftIndex + 1;

            addPair(pairIndex);

            document.getElementById("thumb-" + index).src = images[index];

            updateMainView();
        }

        function connect() {
            let ws = new WebSocket("ws://localhost:"""
        + str(port)
        + """/images");

            ws.onmessage = function(message) {
                let msg = JSON.parse(message.data);
                addImage(msg.index, msg.image);
            };

            ws.onclose = function() {
                setTimeout(connect, 10000);
            };
        }

        connect();
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
        log.info(f"Got connection from {self.request.remote_ip}")
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
            log.info(f"Client {self.request.remote_ip} disconnected")


class HTTPHandler(RequestHandler):
    def get(self):
        # Answer the HTTP request with the index.html page
        self.write(open("index.html").read())


class WebServer:
    def __init__(self, extension, port=8080):
        self.extension = extension
        self.port = port

        # Load the current images
        self.images = {}
        for index in range(MAX_FEEDS * 2):
            self._load_image(f"latest_{index}.{extension}")

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
            raise OSError(f"Port {self.port} is already being used, choose a different one!")

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
        except OSError:
            image_string = ""
        self.images[index] = image_string
        return index, image_string

    def update(self, image_file):
        # Update the image in memory before sending it
        index, image_string = self._load_image(image_file)
        clients = self.connected_clients
        for client in clients:
            self.ioloop.add_callback(client.send_image, index)
