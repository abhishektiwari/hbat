#!/usr/bin/env python
"""
HBAT Web Server launcher script.

This script launches the HBAT web interface.
Usage: ./hbat-web
"""

if __name__ == "__main__":
    from hbat.server.app import create_app

    create_app()
