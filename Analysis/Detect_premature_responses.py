# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 12:34:26 2023

We might also want to get the following info for analyses, to avoid needing to compare this output with the 
behavioral data files:
  - trial (counter for trial)
  - block (based on counter for trial)
  - FP (can be computed from time differences)
  - fp n-1 (previous computation stored)
"""

# Import libraries
import pandas as pd
import glob
import os


# Instruction texts for comparisons
action_start = 'instrucoes_practica_action: text = Nesta primeira parte da tarefa, a cruz preta aparecerá na tela e permanecerá lá até você pressionar a barra de espaço para fazê-la mudar de cor. Utilize sua mão não-dominante para pressionar a barra de espaço (p. ex., se você é destro, utilize a mão esquerda para pressionar a barra de espaço). Só então a cruz mudará de preta para verde. Nesse momento, você deverá se preparar para a apresentação do estímulo.'
external_start = 'instrucoes_pratica_external: text = Nesta primeira parte da tarefa, a cruz mudará de cor automaticamente. Você deverá ficar atento ao momento em que a mudança de cor acontecer e se preparar para a apresentação do estímulo.'

action_middle = 'instrucoes_practica_action: text = Nos próximos blocos, a tarefa será um pouco diferente. Ao invés de esperar a cruz mudar de cor sozinha, você deverá pressionar a barra de espaço para fazê-la mudar de cor (de preta para verde). Utilize sua mão não-dominante para pressionar a barra de espaço (p. ex., se você é destro, utilize a mão esquerda para pressionar a barra de espaço).'
external_middle = 'instrucoes_pratica_external: text = Nos próximos blocos, a tarefa será um pouco diferente. Ao invés de pressionar a barra de espaço para fazer a cruz mudar de cor, ela mudará de cor sozinha. Você deverá ficar atento ao momento em que a mudança de cor acontecer.'

practice_start = 'Antes de iniciarmos, vamos praticar um pouco. Pressione a barra de espaço para iniciar a fase de prática.'
practice_middle = 'Primeiro, pratique algumas vezes para ter certeza de que compreendeu a mudança. Pressione a barra de espaço para iniciar a fase de prática.'

# Path to log files
filesPath = '.\Data\Log_files'

# Get list of log files
FileList=glob.glob(filesPath + '/*.log')
FileList.sort()

nFiles=int(len(FileList))

# Dataset to save preResp info
preRespData = pd.DataFrame(columns = ['ID', 'preRespAction', 'preRespExternal', 'preRespTotal'])

for iFile, fileName in enumerate(FileList):
    
    # read file in correct format 
    logData = open(fileName).read()
    
    # Split to read line by line
    lines = logData.split('\n')
    numberLines = len(lines)
    
    # Get participant ID from file name
    ID = fileName[17:20]
    
    # Create placeholder condition variable
    condition = None
    
    # Set practice to False (changes program behavior when set to True)
    practice = False
    
    # premature response counters by condition
    preRespCountAction = 0
    preRespCountExternal = 0
    
    # Read line by line
    for line in range(numberLines - 1):
      # Split line into columns (only the third columns is of interest)
      columns = lines[line].split('\t')
      
      # Perform comparisons on third column
      if practice == False:
        if len(columns) == 1:
          # Check for practice start and set practice to True if started; does nothing otherwise
          if columns[0] in [practice_start, practice_middle]:
            practice = True # does not look for responses
            
        # run through lines with ['Keydown: ArrowRight', 'Keydown: ArrowLeft']
        if len(columns) == 3:
          if columns[2] in ['Keydown: ArrowRight', 'Keydown: ArrowLeft']:
            if condition == 'action':
              if prevLine[2] in ['Fixation_action: autoDraw = true', 'Keydown: \n', 'Fixation_action: autoDraw = null', 'target: image = {}', 'Fixation_2: autoDraw = true', 'Fixation_2: autoDraw = null']:
                preRespCountAction += 1
            if condition == 'external':
              if prevLine[2] in ['Fixation_ext: autoDraw = true', 'Fixation_ext: autoDraw = null', 'target: image = {}', 'Fixation_2: autoDraw = true', 'Fixation_2: autoDraw = null']:
                preRespCountExternal += 1
               
        
      # Check if practice has ended and determine condition
      if practice == True:
        if len(columns) == 1:
          # Check that practice was not repeated
          if columns[0] == 'A fase de prática terminou. Agora iremos iniciar a tarefa. Ela será igual à prática, porém mais longa.':
            practice = False
            
        if len(columns) == 3:
          # Determine condition based on text
          if columns[2] in [external_start, external_middle]:
            condition = 'external'
          elif columns[2] in [action_start, action_middle]:
            condition = 'action'
    
      # Store this line for comparison with the next
      prevLine = columns
    
    # Add to dataset
    partPreData = [ID, preRespCountAction, preRespCountExternal, preRespCountAction + preRespCountExternal]
    preRespData.loc[len(preRespData.index)] = partPreData

preRespData.to_csv('./Analysis/premature_responses.csv')
