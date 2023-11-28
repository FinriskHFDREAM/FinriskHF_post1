## NOTES ON SUBMISSION RETRIVAL FROM SYNAPSE

Details informations on running a challenge hosted in synapse is available here: https://help.synapse.org/docs/Running-a-Challenge.2163409505.html#:~:text=Running%20a%20challenge%20on%20Synapse,the%20infrastructure%20for%20submission%20evaluation

### Download

To list and download submission, please execute this comment.
Please note that the credential to access the Synapse must be acquired.
For details use of accessing synapse using API and setting is available here: https://help.synapse.org/docs/Installing-Synapse-API-Clients.1985249668.html

python ./challenge_infrastructure/SynListDown_final.py -d <directory of the submission>

Using the script, the system will download submission that is not available yet in the target system from synapse. The script also will automatically update the status of submissions.



### Runing

After downloading the submission, all model is subsequently evaluated using the real dataset according to the participants submission.

in case of contanerized image:
```
singularity run your_singularity.sif  <Argument for input path>
```

in case of zipped file with separate singularity image and code, the run instructions provided in `Run` file will be use.

Please see information here for detail: https://www.synapse.org/#!Synapse:syn27130803/wiki/619286

### Scoring

After the model is succesfully run, we evaluated the model to calculate the Harrel's C and hosmer-lemeshow p -value using the command bellow:

Rscript challenge_infrastructure/scoring_evaluation_withrannge_final.R

### Update the metric in leaderboard

Upon completing the evaluations, the status of Leaderboard is updated using command bellow:

python ./challenge_infrastructure/update_status.py -s <submission ID> -d <directory of submission> -n <name of submission>



