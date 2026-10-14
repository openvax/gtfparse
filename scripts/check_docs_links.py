"""Check rendered internal links and fragment targets."""
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit
import yaml


class Page(HTMLParser):
    def __init__(self, text):
        super().__init__()
        self.ids = set()
        self.links = []
        self.feed(text)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if attrs.get("id"):
            self.ids.add(attrs["id"])
        if tag == "a" and attrs.get("href"):
            self.links.append(attrs["href"])


def main():
    root = Path("site").resolve()
    base = urlsplit(yaml.safe_load(Path("mkdocs.yml").read_text())["site_url"]).path.rstrip("/")
    pages = {p.resolve(): Page(p.read_text()) for p in root.rglob("*.html")}
    failures = []
    count = 0
    for source, page in pages.items():
        for href in page.links:
            url = urlsplit(href)
            if url.scheme or url.netloc:
                continue
            path = unquote(url.path)
            if path.startswith("/"):
                if not (path == base or path.startswith(base + "/")):
                    continue
                target = (root / path[len(base):].lstrip("/")).resolve()
            else:
                target = (source.parent / path).resolve() if path else source
            if target.is_dir():
                target /= "index.html"
            count += 1
            if not target.is_relative_to(root) or not target.exists():
                failures.append(f"{source.relative_to(root)}: missing {href}")
            elif url.fragment and target in pages and unquote(url.fragment) not in pages[target].ids:
                failures.append(f"{source.relative_to(root)}: missing fragment {href}")
    if failures:
        raise SystemExit("\n".join(failures))
    print(f"Checked {len(pages)} pages and {count} internal links.")


if __name__ == "__main__":
    main()
