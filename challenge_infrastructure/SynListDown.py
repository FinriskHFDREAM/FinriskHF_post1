#!/usr/bin/env python
#adapted from: https://github.com/Sage-Bionetworks/SynapseWorkflowExample/blob/master/downloadSubmissionFile.cwl
#list submission, download, and create directory accordingly
import synapseclient
syn = synapseclient.Synapse()
syn.login()
import os
from challengeutils import utils
import time
timestr = time.strftime("%Y%m%d")
import pandas as pd
import argparse
from zipfile import ZipFile
parser = argparse.ArgumentParser()
#parser.add_argument("-s", "--submissionId", required=True, help="Submission ID")
#parser.add_argument("-c", "--synapseConfig", required=True, help="credentials file")
parser.add_argument("-d", "--workingdir", required=True, help="working directory path")
args = parser.parse_args()
submissions = syn.tableQuery("SELECT * FROM syn38128821 WHERE status NOT IN ('INVALID', 'CLOSED','SCORED')").asDataFrame()
# Find duplicate submissions based on entityid.
duplicate_subs = submissions.loc[
  submissions
    .sort_values(['entityid', 'createdOn'])  # sort by entity, then by date
    .duplicated(subset=['entityid'])         # get duplicates based on entity value
, 'id'].tolist()  # only keep the `id` values and put them into a list

# Update the status of duplicate submissions to CLOSED.
for sub_id in duplicate_subs:
    utils.change_submission_status(syn, sub_id, "CLOSED")

#listing the ID with RECEIVED status and download the files to the temporary location
#based on submission ID, change the temporary dir locations and ID later
#ID=syn.getSubmissions(9615015,status="RECEIVED")
submissions2 = syn.tableQuery("SELECT * FROM syn38128821 WHERE status NOT IN ('INVALID', 'CLOSED','SCORED')").asDataFrame()
ID = submissions2['id'].tolist()
appended_data = []
print(ID)
for submissionId in ID:
    print(submissionId)
    directory = (args.workingdir+"/"+"submission"+str(submissionId))
    os.makedirs(directory,exist_ok=True)
    sub =  syn.getSubmission(submissionId,downloadLocation=directory)
    file = os.listdir(directory)
    prefix = file[0].split(".")[0]
    suffix = file[0].split(".")[1]
    dirnew = args.workingdir+"/"+prefix
  #print(dirnew)
    if suffix == "sif":
        os.rename(directory,dirnew)
        command="Singularity run"+" "+dirnew+"/"+file[0]
        print(command)
    elif suffix=="zip":
        print(submissionId)
        os.rename(directory,dirnew)
        with ZipFile((dirnew+"/"+file[0]), 'r') as zipObj:
            zipObj.extractall(dirnew)
        content=os.listdir(dirnew)
        print(content)
        command=""
        if "Run" in content:
            with open(dirnew+"/"+"Run") as f:
                lines = f.read().splitlines()
                command=""
                #command="Singularity exec"+" "+lines[2]+" "+ dirnew+"/"+lines[3] + " "+dirnew+"/"+lines[4]
        else:
            #command="command file not found in:"+str(submissionId)+";"+prefix
            command=""
            print("command file not found in:"+str(submissionId)+";"+prefix)
        
    else:
        utils.change_submission_status(syn, submissionId, "INVALID")
        print('Expected ".zip" or ".sif" file,but found other type in:'+str(submissionId)+";"+prefix)
    result = {'entityId':sub.entity.id,'entityVersion':sub.entity.versionNumber,'submission.name':prefix,"submission.id":submissionId,'command':command}
#append to dataframe, save to file with date
    df_result = pd.DataFrame(result,index=[0])
    appended_data.append(df_result)
appended_data = pd.concat(appended_data)
appended_data.to_csv(args.workingdir+"/"+timestr+"_record_Submission"+".csv")

