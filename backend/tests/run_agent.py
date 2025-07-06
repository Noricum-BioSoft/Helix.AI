import asyncio

from backend.agent import handle_command


def main():
    command = "align the following sequences."
    result = asyncio.run(handle_command(command))
    print(result)


if __name__ == "__main__":
    main()
