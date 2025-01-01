import os
from pathlib import Path

from dotenv import load_dotenv


def load_environment():
    """Load environment variables from .env file"""
    env_path = Path(__file__).resolve().parents[2] / ".env"
    if env_path.exists():
        load_dotenv(dotenv_path=env_path)
    else:
        raise FileNotFoundError(".env file not found")

    # Verify critical environment variables are set
    required_vars = ["USERNAME", "PROJECT_CODE", "PROJECT_DIR"]

    missing_vars = [var for var in required_vars if not os.getenv(var)]

    if missing_vars:
        raise EnvironmentError(
            f"Missing required environment variables: {', '.join(missing_vars)}"
        )
