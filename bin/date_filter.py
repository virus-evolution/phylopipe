#!/usr/bin/env python3

import argparse
import datetime
from datetime import date
import csv

def parse_args():
    parser = argparse.ArgumentParser(description="""Filter UK sequences based on metadata""",
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--in-metadata', dest = 'in_metadata', required=True, help='CSV to process')
    parser.add_argument('--out-metadata', dest = 'out_metadata', required=True, help='CSV to write out')
    parser.add_argument('--date', required=False, help='Date to measure time window from, if none provided uses todays date')
    parser.add_argument('--time-window', required=False, default=30, type=int, help='Number of days')
    parser.add_argument('--filter-column', dest = 'filter_column', required=False, default="date_filter", help='Column to update with filter output')
    parser.add_argument('--restrict', action="store_true", help='Only output recent metadata rows')

    args = parser.parse_args()
    return args

def filter_by_date(in_metadata, out_metadata, todays_date, time_window, filter_column, restrict):
    """
    input is CSV, last column being the representative outgroups:
    """
    window = datetime.timedelta(time_window)
    if not todays_date:
        todays_date = date.today()
    todays_date = datetime.datetime.strptime(todays_date, '%Y-%m-%d').date()
    print(window, todays_date)

    with open(in_metadata, 'r', newline = '') as csv_in, \
         open(out_metadata, 'w', newline = '') as csv_out:

        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        fieldnames = reader.fieldnames
        for column in ["why_excluded", filter_column]:
            if column not in fieldnames:
                fieldnames.append(column)
        writer = csv.DictWriter(csv_out, fieldnames = fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()

        for row in reader:
            for column in ["why_excluded", filter_column]:
                if column not in row:
                    row[column] = ""

            if row["why_excluded"] not in [None,"None",""]:
                if not restrict:
                    writer.writerow(row)
                continue

            try:
                date = datetime.datetime.strptime(row["published_date"], '%Y-%m-%d').date()
            except:
                try:
                    date = datetime.datetime.strptime(row["sample_date"], '%Y-%m-%d').date()
                except:
                    row["why_excluded"] = "no sample_date"
                    if not restrict:
                        writer.writerow(row)
                    continue

            if (todays_date - window) > date:
                row[filter_column] = "sample_date older than %s days" %time_window
                if not restrict:
                    writer.writerow(row)
                continue

            writer.writerow(row)


def main():
    args = parse_args()
    filter_by_date(args.in_metadata, args.out_metadata, args.date, args.time_window, args.filter_column, args.restrict)

if __name__ == '__main__':
    main()
