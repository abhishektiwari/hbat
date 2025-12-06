"""
Entry point for HBAT web server.

Run with: python -m hbat.server
"""

from .app import create_app

if __name__ == "__main__":
    create_app()
