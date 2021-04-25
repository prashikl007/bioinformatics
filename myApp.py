# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 23:21:00 2021

@author: PRASHIK
"""

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from PyQt5.QtWidgets import (QApplication, QWidget, QFileDialog, QPushButton, QLabel, QGridLayout, QVBoxLayout)

from Bio import SeqIO
from collections import Counter
from Bio.SeqUtils import molecular_weight
from Bio.SeqUtils import GC
from Bio.SeqUtils import MeltingTemp as mt

import pandas as pd


class MainWindow(QWidget):
    
    def __init__(self):
        super().__init__()
        
        self.setWindowTitle("DNA Sequence Analysis - Prashik Lokhande")
        self.setLayout(QVBoxLayout())
        
        my_label = QLabel("DNA Sequence Analysis from the FASTA Database, (FASTA databse can be found on NCBI website). Build by Prashik Lokhande")
        self.layout().addWidget(my_label)
        
        
        self.visualize()
        self.show()
        
        
    def visualize(self):
        container = QWidget()
        container.setLayout(QGridLayout())
        
        label_1 = QLabel("PLease Select FASTA file")
        
        button_1 = QPushButton("Select file", clicked = lambda: self.get_plot_and_details())
        
        self.record_description = QLabel("Record Description")
        
        gc_content_label = QLabel("GC Content percentage =")
        self.gc_content_field = QLabel("0")
        
        AT_content_label = QLabel("AT Contents percentage =")
        self.AT_content_field = QLabel("0")
        
        mol_weight = QLabel("Molecular weight =")
        self.mol_weight_field = QLabel("0")
        
        melting_point = QLabel("Melting point in celsius =")
        self.melting_point_field = QLabel("0")
        
        longest_acid =QLabel("Longest Amino acid is")
        self.longest_acid = QLabel("Acid")
        
        self.canvas = FigureCanvas(plt.Figure(figsize=(10, 4)))
        self.ax = self.canvas.figure.subplots()
        
        container.layout().addWidget(label_1, 0,0)
        container.layout().addWidget(button_1, 1,0)
        container.layout().addWidget(self.record_description, 2, 0, 1, 1)
        container.layout().addWidget(gc_content_label, 3, 1)
        container.layout().addWidget(self.gc_content_field, 4, 1)
        container.layout().addWidget(AT_content_label, 5, 1)
        container.layout().addWidget(self.AT_content_field, 6, 1)
        container.layout().addWidget(mol_weight, 7, 1)
        container.layout().addWidget(self.mol_weight_field, 8, 1)
        container.layout().addWidget(melting_point, 9, 1)
        container.layout().addWidget(self.melting_point_field, 10, 1)
        container.layout().addWidget(longest_acid, 11, 1)
        container.layout().addWidget(self.longest_acid, 12, 1)
        
        container.layout().addWidget(self.canvas, 3, 0, 15, 1)
        
        self.layout().addWidget(container)
        
    
    def get_plot_and_details(self):
        filepath, _ = QFileDialog.getOpenFileName(self, 'select FASTA file')
        record = SeqIO.read(filepath,"fasta")
        
        self.record_description.setText(str(record.description)) 
        
        dna = record.seq
        mrna = dna.transcribe()
        protein = mrna.translate()
        
        gc = GC(dna)
        self.gc_content_field.setText(str(gc))
        
        AT_content = float(dna.count('A') + dna.count('T'))/len(dna)*100
        self.AT_content_field.setText(str(AT_content))
        
        mol_weight = molecular_weight(dna)
        self.mol_weight_field.setText(str(mol_weight))
        
        melting_point = mt.Tm_GC(dna, strict=False)
        self.melting_point_field.setText(str(melting_point))
        
        pr_freq = Counter(protein)
        self.ax.clear()
        self.ax.bar(pr_freq.keys(), pr_freq.values())
        self.canvas.draw_idle()
        self.ax.set_title("Amino Acid Contents in the sequence (X-axis Amino acids, Y-axis frequency)")
        
        
        clean = protein.split("*")
        clean = [str(i) for i in clean]
        
        df =pd.DataFrame({"amino_acid":clean})
        df['count'] = df['amino_acid'].str.len()
        largest_10 = df.nlargest(10,'count')
        #print("largest Amino Acid with count are ", largest_10)
        self.longest_acid.setText(str(largest_10))
        

app = QApplication([])
mw = MainWindow()

app.exec_()
