def read_fastq_record(file):
    header = file.readline()

    if not header:
        return None

    sequence = file.readline()
    plus = file.readline()
    quality = file.readline()

    return header, sequence, plus, quality


def write_fastq_record(file, record):
    file.writelines(record)


def gc_content(sequence):
    sequence = sequence.strip().upper()

    if not sequence:
        return 0

    gc_count = sequence.count("G") + sequence.count("C")

    return gc_count / len(sequence) * 100


def average_quality(quality):
    quality = quality.strip()

    if not quality:
        return 0

    quality_scores = [ord(char) - 33 for char in quality]

    return sum(quality_scores) / len(quality_scores)


def check_gc(sequence, gc_bounds):
    gc = gc_content(sequence)

    if isinstance(gc_bounds, (int, float)):
        return gc < gc_bounds

    lower, upper = gc_bounds

    return lower <= gc <= upper


def check_length(sequence, length_bounds):
    length = len(sequence.strip())

    if isinstance(length_bounds, (int, float)):
        return length < length_bounds

    lower, upper = length_bounds

    return lower <= length <= upper


def check_quality(quality, quality_threshold):
    return average_quality(quality) >= quality_threshold


def validate_bounds(bounds):
    if isinstance(bounds, (int, float)):
        return

    if not isinstance(bounds, tuple):
        raise TypeError("Границы должны быть числом или кортежем")

    if len(bounds) != 2:
        raise ValueError("Кортеж границ должен содержать два значения")

    lower, upper = bounds

    if lower > upper:
        raise ValueError(
            "Нижняя граница не может быть больше верхней границы"
        )


def validate_quality_threshold(threshold):
    if not isinstance(threshold, (int, float)):
        raise TypeError("Порог качества должен быть числом")

    if threshold < 0:
        raise ValueError("Порог качества не может быть отрицательным")


