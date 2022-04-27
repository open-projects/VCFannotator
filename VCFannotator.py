#!/usr/bin/env python3

import re
import os
import sys
import argparse
import sqlite3
import random
import string
import gzip


class VCFdb:
    _drop_table = " DROP TABLE IF EXISTS vcf; "
    _create_table = """ CREATE TABLE vcf (
                                    n_rec INT NOT NULL,
                                    n_alt INT NOT NULL,
                                    chrom TEXT,
                                    pos INT NOT NULL,
                                    snp_id TEXT,
                                    ref TEXT,
                                    alt TEXT,
                                    qual TEXT,
                                    filter TEXT,
                                    info TEXT,
                                    format TEXT,
                                    samples TEXT
                                ); """
    _create_index_1 = " CREATE INDEX chr_pos ON vcf (pos, chrom); "
    _create_index_2 = " CREATE INDEX n_rec ON vcf (n_rec); "
    _insert = """ INSERT INTO vcf (n_rec, n_alt, chrom, pos, snp_id, ref, alt, qual, filter, info, format, samples)
                            VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?); """

    def __init__(self, file_name=''):
        self._tmp_db_file = ''
        self._db_file = file_name
        self._connect, self._cursor = self._initialize_db()

    def __del__(self):
        try:
            self._connect.close()
            if os.path.exists(self._tmp_db_file):
                os.remove(self._tmp_db_file)
        except sqlite3.Error as er:
            print(er, file=sys.stderr)
            print("ERROR - Can't close db connection for {} file".format(self._tmp_db_file))

    def _initialize_db(self):
        if len(self._db_file):
            try:
                if os.path.exists(self._db_file):
                    connect = sqlite3.connect(self._db_file)
                    cursor = connect.cursor()

                    return connect, cursor  # it will use earlier created db
                else:
                    connect = sqlite3.connect(self._db_file)
                    cursor = connect.cursor()
            except sqlite3.Error as er:
                print(er, file=sys.stderr)
                exit("Can't connect to db {}\n".format(self._db_file))
        else:
            self._tmp_db_file = ''.join(random.choice(string.ascii_lowercase) for i in range(8)) + '.db'
            try:
                connect = sqlite3.connect(self._tmp_db_file)
                cursor = connect.cursor()
            except sqlite3.Error as er:
                print(er, file=sys.stderr)
                exit("Can't create file and connect to db {}\n".format(self._tmp_db_file))

        # cursor.execute(self._drop_table)
        cursor.execute(self._create_table)
        connect.commit()

        return connect, cursor

    def add(self, data_array):  # [n_rec, n_alt, chrom, pos, snp_id, ref, alt, qual, filter, info, format, samples]
        try:
            self._connect.executemany(self._insert, data_array)
            self._connect.commit()
        except sqlite3.Error as er:
            print(er, file=sys.stderr)
            exit("ERROR - Can't commit to db {}\n".format(' '.join(str(e) for e in data_array)))

        return self._cursor.lastrowid

    def index(self):
        self._cursor.execute(self._create_index_1)
        self._cursor.execute(self._create_index_2)
        self._connect.commit()

    def attach(self, db_name, db_file):
        try:
            self._cursor.execute(" ATTACH DATABASE ? AS ?; ", (db_file, db_name))
        except sqlite3.Error as er:
            print(er, file=sys.stderr)
            exit("ERROR - Can't attach the database file {}\n".format(db_file))

        return self._cursor

# end of 'VCFdb' class


class VCFrec:
    _freq_pattern = r'FREQ=([^;]+)'

    def __init__(self, record):
        self.chrom = '.'
        self.pos = 0
        self.snp_id = '.'
        self.ref = '.'
        self.alt = ()
        self.qual = '.'
        self.filter = '.'
        self.info = '.'
        self.format = ''
        self.samples = ()

        if re.match(r'^#', record):
            raise ValueError("ERROR - Commented line: '{}'".format(record))
        else:
            fields = re.split(r'\t', record.strip())
            if len(fields) > 7:
                self.chrom, self.pos, self.snp_id, self.ref, alt, self.qual, self.filter, self.info = fields[0:8]
                alt = re.sub(r'\*', '.', alt)

                if re.search(r'[^ATGCatgc]', self.ref):
                    print("WARNING - Wrong symbol in reference sequence: {}".format(self.ref))
                if re.search(r'[^ATGCatgc,.]', alt):
                    print("WARNING - Wrong symbol in alternative sequence: {}".format(alt))

                self.alt = re.split(r',', alt)
                if len(fields) > 9:
                    self.format = fields[8]
                    self.samples = fields[9:]
            else:
                raise ValueError("ERROR - Wrong format of the input record: '{}'".format(record))

    def get_string(self):
        record = '\t'.join((
            self.chrom,
            str(self.pos),
            self.snp_id,
            self.ref,
            ','.join(self.alt),
            self.qual,
            self.filter,
            self.info
        ))

        if self.format:
            record += '\t' + self.format
            if len(self.samples):
                record += '\t' + '\t'.join(self.samples)
            else:
                raise ValueError("ERROR - No sample data: '{}'".format(record))

        return record

    def mod_info_freq(self, n_alt):
        re_freq = re.search(r'FREQ=([^;]+)', self.info)
        mod_freq = ''
        if re_freq:
            freq_group = re_freq.group(1)
            for pop_freq in freq_group.split('|'):
                freq_ver, freq_set = re.split(r':', pop_freq)
                freq_array = re.split(r',', freq_set)
                if len(freq_array) < 2 or len(freq_array) != len(self.alt) + 1:
                    raise ValueError("ERROR - Wrong allele frequency data: '{}'".format(self.info))

                ref_freq = freq_array[0]
                alt_freq = freq_array[n_alt]

                if mod_freq:
                    mod_freq += '|' + freq_ver + ':' + ref_freq + ',' + alt_freq
                else:
                    mod_freq = 'FREQ=' + freq_ver + ':' + ref_freq + ',' + alt_freq

            mod_info = re.sub(r'FREQ=[^;]+', mod_freq, self.info)

            return mod_info

        return self.info

    def get_array(self):
        array = []
        n_alt = 1
        if len(self.alt) > 1:
            for alt in self.alt:
                array.append([
                    n_alt,
                    self.chrom,
                    self.pos,
                    self.snp_id,
                    self.ref,
                    alt,
                    self.qual,
                    self.filter,
                    self.mod_info_freq(n_alt),
                    self.format,
                    '\t'.join(self.samples)
                ])
                n_alt += 1
        else:
            array.append([
                n_alt,
                self.chrom,
                self.pos,
                self.snp_id,
                self.ref,
                self.alt[0],
                self.qual,
                self.filter,
                self.info,
                self.format,
                '\t'.join(self.samples)
            ])

        return array

    def add_info(self, inf_key, inf_value):
        if inf_key and inf_value:
            return self.info.rstrip(';') + ';{}={}'.format(inf_key, inf_value)
        raise ValueError("ERROR - Bad info values: {}={}".format(inf_key, inf_value))

# end of 'VCFrec' class


class VCF:
    def __init__(self, vcf_file_name):
        self._vcf_file = vcf_file_name
        self._header, self._columns = self._header()
        self._insert_chunk_size = 100000
        self._db = None

    def _header(self):
        hdr, cln = '', ''
        with gzip.open(self._vcf_file, 'rt') as f:
            for line in f:
                if re.match(r'##', line):
                    hdr += line
                elif re.match(r'#CHROM', line):
                        cln = line
                        break
                elif len(hdr) > 0:
                    break
        if len(hdr) == 0 or len(cln) == 0:
            exit("ERROR: bad VCF header in file: {}\n".format(self._vcf_file))

        return hdr, cln

    def load2db(self, db_file_name=''):
        db = VCFdb(db_file_name)

        n_rec = 0
        data = []
        with gzip.open(self._vcf_file, 'rt') as f:
            for line in f:
                if len(data) > self._insert_chunk_size:
                    db.add(data)
                    data = []
                line = line.strip()
                if re.match(r'^[^#]', line):
                    record = VCFrec(line)
                    n_rec += 1
                    for data_array in record.get_array():
                        data.append([n_rec] + data_array)
            if len(data):  # the last chunk
                db.add(data)

        db.index()

        self._db = db
        return n_rec

    def get_header(self):
        return self._header

    def get_db(self):
        return self._db

    def annotation2csv(self, an_db_file, out_file='annotation.csv'):
        select = """ SELECT t1.*, t2.snp_id, t2.ref, t2.alt, t2.info FROM vcf AS t1 
                    LEFT JOIN annotation.vcf AS t2 ON t1.pos = t2.pos AND t1.chrom = t2.chrom AND
                    (t1.ref = t2.ref AND t1.alt = t2.alt OR t1.ref = t2.alt AND t1.alt = t2.ref)
                    ORDER BY t1.n_rec, t1.n_alt; """

        header = ('#_rec', 'n_alt', 'chrom', 'pos', 'var_id', 'ref', 'alt', 'qual', 'filter', 'info', 'format',
                  'samples', 'snp_id', 'snp_ref', 'snp_alt', 'snp_info')

        n = 0
        try:
            cursor = self._db.attach('annotation', an_db_file)
            with open(out_file, 'w') as csv:
                csv.write('\t'.join(header) + "\n")
                for row in cursor.execute(select):
                    csv.write('\t'.join(re.sub(r'^None$', r'\\N', str(field)) for field in row) + "\n")
                    n += 1
        except sqlite3.Error as er:
            print(er, file=sys.stderr)
            exit("ERROR - Can't make the annotation... 8(\n")

        return n

# end of 'VCF' class


def main():
    parser = argparse.ArgumentParser(description="A program to annotate VCF files using VCF files.")
    parser.add_argument('-i', help='input - gzipped VCF file to annotate')
    parser.add_argument('-a', help='annotation - gzipped VCF file with annotation (TOPMED, gnomAD, etc.)')
    parser.add_argument('-o', default='output', help='output file')

    args = parser.parse_args()
    in_file = args.i
    an_file = args.a
    out_file = args.o

    an_db_file = an_file + '.db'
    if not os.path.exists(an_db_file):
        an_vcf = VCF(an_file)
        an_vcf.load2db(an_db_file)
        an_vcf = None  # to destroy the internal database object and release the database file

    in_vcf = VCF(in_file)
    in_vcf.load2db()

    n_annotated_rows = in_vcf.annotation2csv(an_db_file, out_file)

    print('{} rows in the output file\ndone...'.format(n_annotated_rows))


if __name__ == "__main__":
    main()
