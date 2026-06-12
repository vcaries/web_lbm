# Dev server with real COOP/COEP headers (crossOriginIsolated without the
# service-worker shim, so the multithreaded engine starts with no reload).
# Run from the repo root:  python serve_coi.py  ->  http://localhost:8001/app/
import http.server


class Handler(http.server.SimpleHTTPRequestHandler):
    extensions_map = {
        **http.server.SimpleHTTPRequestHandler.extensions_map,
        ".js": "text/javascript",
        ".mjs": "text/javascript",
        ".wasm": "application/wasm",
    }

    def end_headers(self):
        self.send_header("Cross-Origin-Opener-Policy", "same-origin")
        self.send_header("Cross-Origin-Embedder-Policy", "require-corp")
        self.send_header("Cross-Origin-Resource-Policy", "same-origin")
        super().end_headers()


if __name__ == "__main__":
    http.server.test(HandlerClass=Handler, port=8001)
