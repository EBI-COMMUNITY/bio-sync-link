# std
import time
# custom
import ena_gbif_mapper as egm


if __name__ == "__main__":
    # create a df with test data
    df = egm.read_ena_dump_pd('input/ena_dump')
    start = time.time()
    res = egm.process(df)
    end = time.time()
    print("Time elapsed: ", end - start)
