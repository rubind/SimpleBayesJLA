try:
    from commands import getoutput
except:
    from subprocess import getoutput
import time

def count_jobs():
    job_count = int(getoutput("squeue | grep " + whoami + " | wc").split(None)[0])
    print("job_count", job_count)
    return job_count
    
whoami = getoutput("whoami")

while count_jobs() > 1:
    print("Waiting...")
    time.sleep(900)
