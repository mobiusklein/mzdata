import argparse
import base64
from urllib import parse as urlparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("img_path")
    parser.add_argument("tag_name")
    return parser.parse_args()


def make_data_uri_from_path(path):
    content = open(path, "rb").read()
    content = base64.b64encode(content)
    content = urlparse.quote(content)
    return f"data:image/png;base64,{content}"


def main():
    args = parse_args()
    uri = make_data_uri_from_path(args.img_path)
    print(f" [{args.tag_name}]: {uri}")


if __name__ == "__main__":
    main()
