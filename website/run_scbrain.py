"""scBRAIN local launcher: starts the site on this computer and opens it in the browser.
Data are read from the 'data' folder next to this program (or next to run_scbrain.py when run with Python).
Close the window to stop."""
import os, socket, sys, threading, time, webbrowser

BASE = os.path.dirname(sys.executable if getattr(sys, 'frozen', False) else os.path.abspath(__file__))
os.environ.setdefault('SCBRAIN_DATA', os.path.join(BASE, 'data').replace('\\', '/'))
if getattr(sys, 'frozen', False):
    sys.path.insert(0, sys._MEIPASS)
else:
    sys.path.insert(0, BASE)


def free_port(start=8765):
    for p in range(start, start + 50):
        with socket.socket() as s:
            try:
                s.bind(('127.0.0.1', p))
                return p
            except OSError:
                continue
    raise SystemExit('No free port found.')


def main():
    if not os.path.isdir(os.path.join(os.environ['SCBRAIN_DATA'], 'tables')):
        print('Data folder not found:', os.environ['SCBRAIN_DATA'])
        input('Press Enter to close.')
        return
    print('scBRAIN is starting, please wait (about 10-30 seconds the first time)...')
    import uvicorn
    from app.main import app
    port = free_port()
    url = f'http://127.0.0.1:{port}'

    def opener():
        for _ in range(120):
            try:
                with socket.create_connection(('127.0.0.1', port), timeout=1):
                    break
            except OSError:
                time.sleep(0.5)
        if not os.environ.get('SCBRAIN_NO_BROWSER'):
            webbrowser.open(url)
        print(f'\nscBRAIN is running at {url}\nIt has been opened in your browser. Keep this window open while you use it; close it to stop.\n')

    threading.Thread(target=opener, daemon=True).start()
    uvicorn.run(app, host='127.0.0.1', port=port, log_level='warning')


if __name__ == '__main__':
    main()
