import plac, sys

from tinrun import tin, overlap, bam2gtf, tin_test, tin_prev, overlap_prev

SUB_COMMANDS = {'tin': tin.run, 'overlap': overlap.run, 'bam2gtf': bam2gtf.run,
                'tin_test': tin_test.run, 'tin_prev': tin_prev.run, 'overlap_prev': overlap_prev.run}

USAGE = f"""
   tincheck: Check the coverage evenness of a transcript \n

   tincheck tin       : calculate TIN
   tincheck overlap   : calculate transcript overlap
   tincheck bam2gtf   : create a GTF file from bam file
   tincheck tin_test  : calculate TIN (for tests)
   tincheck tin_prev  : calculate TIN (previos version, with read length, for tests)
   tincheck overlap_prev : calculate transcript overlap (previous version with read length, for tests)

   Run each command for more help.
   """


def run():
    """Calculate Transcript Integrity Number (TIN)"""

    # Print usage when no parameters are passed.
    if len(sys.argv) == 1:
        print(USAGE)
        sys.exit(1)

    if len(sys.argv) < 2:
        print("Enter a subcommand")
        sys.exit(1)

    cmd = sys.argv[1]

    sys.argv.remove(cmd)

    # Raise an error is not a valid subcommand.
    if cmd not in SUB_COMMANDS:
        print(USAGE, file=sys.stderr)
        print(f"invalid command: {cmd}")
        sys.exit(-1)

    func = SUB_COMMANDS[cmd]
    plac.call(func)


if __name__ == '__main__':
    run()
