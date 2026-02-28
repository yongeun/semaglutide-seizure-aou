#!/usr/bin/env python
# coding: utf-8

# In[1]:


# The codes are meanto run sequentially from [0].py -> [1].py -> [2].py -> [3].py -> [4]. py
# 1. Building Reference Dataframe - Python Code

import os
import pandas as pd
import subprocess
import numpy as np

ehr_query = f"""
with ehr as (
    select
        p.person_id,
        p.birth_datetime as date_of_birth,
        c_race.concept_name as race,
        c_sex.concept_name as sex,
        c_ethn.concept_name as ethnicity
    from `{os.environ['WORKSPACE_CDR']}.person` p
    left join `{os.environ['WORKSPACE_CDR']}.concept` c_race
         on p.race_concept_id = c_race.concept_id
    left join `{os.environ['WORKSPACE_CDR']}.concept` c_sex
         on p.sex_at_birth_concept_id = c_sex.concept_id
    left join `{os.environ['WORKSPACE_CDR']}.concept` c_ethn
         on p.ethnicity_concept_id = c_ethn.concept_id
    where p.person_id in (
        select distinct person_id
        from `{os.environ["WORKSPACE_CDR"]}.cb_search_person`
        where has_ehr_data = 1
    )
)
select
    ehr.person_id,
    ehr.date_of_birth,
    ehr.race,
    ehr.ethnicity,
    ehr.sex,
    ins1.aname as ins1,
    ins2.aname as ins2,
    obs1.aname as aname1,
    obs2.aname as aname2,
    obs3.aname as aname3,
    obs4.aname as aname4,
    obs5.aname as aname5,
    obs6.aname as aname6,
    obs7.aname as aname7,
    obs8.aname as aname8
from ehr
left join (
    select
        o.person_id,
        answer.concept_name as aname
    from `{os.environ['WORKSPACE_CDR']}.observation` o
    left join `{os.environ['WORKSPACE_CDR']}.concept` answer
         on o.value_source_concept_id = answer.concept_id
    where o.observation_source_concept_id = 1585386
) ins1 on ehr.person_id = ins1.person_id
left join (
    select
        o.person_id,
        answer.concept_name as aname
    from `{os.environ['WORKSPACE_CDR']}.observation` o
    left join `{os.environ['WORKSPACE_CDR']}.concept` answer
         on o.value_source_concept_id = answer.concept_id
    where o.observation_source_concept_id = 43528428
) ins2 on ehr.person_id = ins2.person_id
left join (
    select
        o.person_id,
        answer.concept_name as aname
    from `{os.environ['WORKSPACE_CDR']}.observation` o
    left join `{os.environ['WORKSPACE_CDR']}.concept` answer
         on o.value_source_concept_id = answer.concept_id
    where o.observation_source_concept_id = 1585375
) obs1 on ehr.person_id = obs1.person_id
left join (
    select
        o.person_id,
        answer.concept_name as aname
    from `{os.environ['WORKSPACE_CDR']}.observation` o
    left join `{os.environ['WORKSPACE_CDR']}.concept` answer
         on o.value_source_concept_id = answer.concept_id
    where o.observation_source_concept_id = 1585940
) obs2 on ehr.person_id = obs2.person_id
left join (
    select
        o.person_id,
        answer.concept_name as aname
    from `{os.environ['WORKSPACE_CDR']}.observation` o
    left join `{os.environ['WORKSPACE_CDR']}.concept` answer
         on o.value_source_concept_id = answer.concept_id
    where o.observation_source_concept_id = 1585857
) obs3 on ehr.person_id = obs3.person_id
left join (
    select
        o.person_id,
        answer.concept_name as aname
    from `{os.environ['WORKSPACE_CDR']}.observation` o
    left join `{os.environ['WORKSPACE_CDR']}.concept` answer
         on o.value_source_concept_id = answer.concept_id
    where o.observation_source_concept_id = 1585860
) obs4 on ehr.person_id = obs4.person_id
left join (
    select
        o.person_id,
        answer.concept_name as aname
    from `{os.environ['WORKSPACE_CDR']}.observation` o
    left join `{os.environ['WORKSPACE_CDR']}.concept` answer
         on o.value_source_concept_id = answer.concept_id
    where o.observation_source_concept_id = 1586198
) obs5 on ehr.person_id = obs5.person_id
left join (
    select
        o.person_id,
        answer.concept_name as aname
    from `{os.environ['WORKSPACE_CDR']}.observation` o
    left join `{os.environ['WORKSPACE_CDR']}.concept` answer
         on o.value_source_concept_id = answer.concept_id
    where o.observation_source_concept_id = 1586201
) obs6 on ehr.person_id = obs6.person_id
left join (
    select
        o.person_id,
        answer.concept_name as aname
    from `{os.environ['WORKSPACE_CDR']}.observation` o
    left join `{os.environ['WORKSPACE_CDR']}.concept` answer
         on o.value_source_concept_id = answer.concept_id
    where o.observation_source_concept_id = 1586207
) obs7 on ehr.person_id = obs7.person_id
left join (
    select
        o.person_id,
        answer.concept_name as aname
    from `{os.environ['WORKSPACE_CDR']}.observation` o
    left join `{os.environ['WORKSPACE_CDR']}.concept` answer
         on o.value_source_concept_id = answer.concept_id
    where o.observation_source_concept_id = 1586213
) obs8 on ehr.person_id = obs8.person_id
"""

# Load data from BigQuery using the updated query
ehr_df = pd.read_gbq(ehr_query, dialect="standard", use_bqstorage_api=True)

# Display the first few rows of the DataFrame
print(ehr_df.head())
print(f"Total number of rows: {len(ehr_df)}")

# Check for duplicates by person_id (all in lower case)
duplicate_count = len(ehr_df) - len(ehr_df.drop_duplicates(subset=['person_id']))
print(f"Found {duplicate_count} duplicate rows (based on person_id)")

# Remove duplicates to ensure exactly one row per person_id
print(f"Original row count: {len(ehr_df)}")
ehr_df = ehr_df.drop_duplicates(subset=['person_id'])
print(f"After removing duplicates: {len(ehr_df)} rows")

my_dataframe = ehr_df   
destination_filename = 'ehr_df.csv'
my_dataframe.to_csv(destination_filename, index=False)
my_bucket = os.getenv('WORKSPACE_BUCKET')
args = ["gsutil", "cp", f"./{destination_filename}", f"{my_bucket}/data/"]
output = subprocess.run(args, capture_output=True)
output.stderr

