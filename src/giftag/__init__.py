"""giftag: annotate genomes with exactly the markers gifter evaluates."""

__version__ = "0.1.0"

# The layout of a database directory written by `giftag build`. `annotate`
# refuses a directory whose format it does not know.
DB_FORMAT = 1


class GiftagError(Exception):
    """An error the command line reports as a message rather than a traceback."""
