"""
Author : Yaochieh Yao
This code is for assignment 1 of the BINF6400 Genomic course.
The homework requires assembling sequence reads into contigs
and takes an input file that reads and outputs a file with a
FASTA sequence of contigs labeled numerically. Finally, a
visualization of sequencing depth / read density map is
included, showing the coverage of alignment across each contig.
"""
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler


def read_reads(filename):
    """The function reads txt file and return seqs to a list"""
    with open(filename, "r") as f:
        # Initialize an empty list to store the lines
        seq_list = []

        # Iterate through the file line by line
        for line in f:
            # Remove the trailing newline character '\n' using strip() and append the line to the list
            seq_list.append(line.strip())
    return seq_list


def assemble_reads(contig, seq_list, k):
    """ 
    The function uses k parameter to check overlap from ahead and
    behind position with the matched sequences and align with the
    one with maximum lengths ahead or behind from the overlap.
    """
    counts = []
    while True:
        # Use k base pair ahead and behind to search matched sequence in the seq_list
        match_ahead = contig[:k]
        match_behind = contig[-k:]

        matched_ahead = [seq for seq in seq_list if match_ahead in seq]
        matched_behind = [seq for seq in seq_list if match_behind in seq]

        # Pick up one with the maxium letters from match_ahead or match_behind sequence from each list
        max_letters_before_match = max(
            matched_ahead, key=lambda seq: seq.index(match_ahead), default="")
        max_letters_after_match = max(matched_behind, key=lambda seq: len(
            seq) - seq.index(match_behind) - len(match_behind), default="")

        # Then we use index slicing to get the letters ahead and behind
        ahead = max_letters_before_match[:max_letters_before_match.index(
            match_ahead)] if match_ahead in max_letters_before_match else ""
        behind = max_letters_after_match[-(len(max_letters_after_match) - max_letters_after_match.rfind(
            match_behind) - len(match_behind)):] if match_behind in max_letters_after_match else ""

        # Visualization chunk
        # Pick up one with the maxium letters from match_ahead or match_behind sequence from each list
        letter_counts = [len(seq) - seq.index(match_behind) -
                         len(match_behind) for seq in matched_behind]

        # This part handling the sequence aligngment in the ending when no more letters after the overlap
        if letter_counts:
            max_letters_after_match = max(letter_counts)
            max_index = letter_counts.index(max_letters_after_match)
            behind = matched_behind[max_index][-max_letters_after_match:]
        else:
            # Return 0 letter fater match and max_index is None, behind sequence is empty
            max_letters_after_match = 0
            max_index = None
            behind = ""

        # This code appends the index number in global contig range, and use to calculate the depth
        for letter_count, ele in zip(letter_counts, matched_behind):
            start_numbers = [i for i in range(
                len(contig[:2*-k])+letter_count, len(contig[:2*-k])+len(ele)+letter_count)]
            end_numbers = [i for i in range(
                len(contig)-len(ele), len(contig)-len(ele)+len(ele))]
            if letter_count == 0 and max_letters_after_match == 0:
                counts.append(end_numbers)
            else:
                counts.append(start_numbers)

        # Now, we can remove the elements of matched list from the main seq_list
        seq_list = [
            seq for seq in seq_list if seq not in matched_ahead and seq not in matched_behind]

        # If there are letters ahead and behind then we extend the contig and continue the while loop
        # otherwise check the final contg before output
        if ahead == "" and behind == "":
            # Check the tail and remove any duplicate sequence at the end (if any)
            # From my observation, this cause from k parameters and shorter seuqences (L<20)
            for i in range(k+1, 20):
                if contig[-i:] == contig[-2*i:-i]:
                    contig = contig[:-i]
            return contig, seq_list, counts
        else:
            contig = ahead + contig + behind


def write_contig(filename, contigs):
    """
    The function outputs a file with a FASTA sequence of contigs
    labeled numerically
    """
    with open(filename, "w") as w:
        for i in range(len(contigs)):
            w.write(f" >{i+1}\n{contigs[i]}\n")

    return filename


if __name__ == "__main__":
    """Main Business Logic"""
    # Sequence Assembly Section
    k = 10
    seq_list = read_reads("seqReadFile.txt")
    contigs = []
    all_counts = []
    # Assemble reads until seq_list become empty
    while seq_list:
        contig, seq_list, counts = assemble_reads(seq_list[0], seq_list, k)
        contigs.append(contig)
        all_counts.append(counts)

    # Write a Fasta file for all contigs and label numerially
    fasta_file = "yaochieh_yao.fasta"
    write_contig(fasta_file, contigs)

    # Visualization Section
    # Count how many sequence at each position of contig as seq_depth
    seq_depth = []
    for contig, count_list in zip(contigs, all_counts):
        contig_depth = [0] * len(contig)
        for idx_list in count_list:
            for idx in idx_list:
                contig_depth[idx] += 1
        seq_depth.append(contig_depth)

    # Plot each contig separately
    for contig_depth in seq_depth:
        # Create indices for each position in the contig
        idx = list(range(len(contig_depth)))
        plt.plot(idx, contig_depth,
                 label=f'Contig {seq_depth.index(contig_depth) + 1}')

    # Add labels. title and a legend to distinguish contigs
    plt.xlabel('Position')
    plt.ylabel('Depth')
    plt.title('Depth Plot for Each Contig')
    plt.legend()
    plt.show()

    # # Normzlized version of depth data
    # for contig_depth in seq_depth:
    #     # Normalize contig_depth using Min-Max scaling
    #     scaler = MinMaxScaler()
    #     normalized_depth = scaler.fit_transform(
    #         [[depth] for depth in contig_depth])

    #     # Create indices for each position in the contig
    #     idx = list(range(len(contig_depth)))

    #     # Plot the normalized depth
    #     plt.plot(idx, normalized_depth,
    #              label=f'Contig {seq_depth.index(contig_depth) + 1}')

    # plt.title('Normalized Depth Plot for Each Contig')
    # plt.legend()
    # plt.xlabel('Position')
    # plt.ylabel('Normalized Depth')
    # plt.show()
