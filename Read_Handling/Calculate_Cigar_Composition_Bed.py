import argparse
import re
import pysam


def parse_cigar(cigar):
    pattern = re.compile(r'(\d+)([MIDNSHPX=])')
    operations = re.findall(pattern, cigar)

    operation_counts = {}
    total_length = 0

    for count, op in operations:
        count = int(count)
        if op in operation_counts:
            operation_counts[op] += count
        else:
            operation_counts[op] = count

        if op in "MDN=X":
            total_length += count

    return operation_counts, total_length


def calculate_percentages(operation_counts, total_length):
    percentages = {op: (count / total_length) * 100 for op, count in operation_counts.items()}
    return percentages


def determine_match_flag(cigar):
    if re.fullmatch(r'(\d+M)+', cigar):
        return 2  # Gap less match
    elif re.fullmatch(r'(\d+[MID])+', cigar):
        return 1  # Gapped match
    else:
        return 0  # Otherwise


def process_bam_file(input_file, output_file):
    bamfile = pysam.AlignmentFile(input_file, "rb")
    with open(output_file, 'w') as outfile:
        for read in bamfile:
            if read.is_unmapped:
                continue

            read_name = read.query_name
            chromosome = bamfile.get_reference_name(read.reference_id)
            position = read.reference_start + 1  # 1-based position
            cigar = read.cigarstring

            if cigar is None:
                continue

            operation_counts, total_length = parse_cigar(cigar)
            # check if total_length is 1) a number 2) greater than 0, if not, print read name to STDOUT and skip this read
            if not total_length or not isinstance(total_length, int):
                print(f"Read {read_name} has total_length {total_length} and is skipped.")
                continue

            percentages = calculate_percentages(operation_counts, total_length)
            match_flag = determine_match_flag(cigar)

            output_fields = [read_name, chromosome, str(position)]
            for op in "MIDNSHPX=":
                output_fields.append(f"{percentages.get(op, 0):.2f}")
            output_fields.append(str(total_length))
            output_fields.append(str(match_flag))

            outfile.write('\t'.join(output_fields) + '\n')


def main():
    parser = argparse.ArgumentParser(
        description="Calculate CIGAR code statistics from a BAM file.",
        epilog="""Output format:
        <read_name>\t<chromosome>\t<position>\t<percent_M>\t<percent_I>\t<percent_D>\t<percent_N>\t<percent_S>\t<percent_H>\t<percent_P>\t<percent_X>\t<percent_=>\t<total_alignment_length>\t<match_flag>

        Example:
        NB501801:498:HG53JBGXK:1:11304:26729:3030\tchr1\t9968\t95.00\t1.25\t2.50\t0.00\t1.25\t0.00\t0.00\t0.00\t0.00\t80\t2
        """
    )
    parser.add_argument("input", help="Input BAM file")
    parser.add_argument("output", help="Output file")

    args = parser.parse_args()

    process_bam_file(args.input, args.output)


if __name__ == "__main__":
    main()
