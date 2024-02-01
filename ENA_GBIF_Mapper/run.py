# std
import time
# custom
import ena_gbif_mapper as egm


if __name__ == "__main__":
    sample_size: int = int(input("Enter sample size: "))
    start = time.time()
    egm.process('input/ena_dump', sample_size)
    end = time.time()
    print("Time elapsed: ", end - start)
