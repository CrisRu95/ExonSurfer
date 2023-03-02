#!/usr/bin/env python3
# -*- coding: utf-8 -*-

design_dict = {
        
        # primer paramters:
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 27,
        
        'PRIMER_OPT_TM': 60,
        'PRIMER_MIN_TM': 57,
        'PRIMER_MAX_TM': 63,
        'PRIMER_TM_FORMULA': 1, # use santaLucia (1998) formula 
        
        'PRIMER_OPT_GC_PERCENT': 50,
        'PRIMER_MIN_GC': 20,
        'PRIMER_MAX_GC': 80,
        
        'PRIMER_MAX_POLY_X': 5,
        'PRIMER_GC_CLAMP': 1, 
        'PRIMER_MAX_END_GC': 4,

        # not shown to user
        'PRIMER_WT_GC_PERCENT_GT':1,
        'PRIMER_WT_GC_PERCENT_LT': 1,
                
        # product parameters
        'PRIMER_PRODUCT_OPT_TM': 80,
        'PRIMER_PRODUCT_MIN_TM': 76,
        'PRIMER_PRODUCT_MAX_TM': 90,
        
        'PRIMER_PRODUCT_SIZE_RANGE': [[120, 250]],
        
        # PCR parametres
        'PRIMER_SALT_DIVALENT': 1.5,
        'PRIMER_SALT_MONOVALENT': 50,
        'PRIMER_SALT_CORRECTIONS': 1, # use santaLucia (1998) formula 
        'PRIMER_DNTP_CONC': 0.6, 
        } 