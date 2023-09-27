##### See: https://docs.crc.nd.edu/new_user/quick_start.html

"""
Write default job script settings to a given file "f".
"""
function write_job_array_settings_crc(f; n::Int64=1)
    println(f, "#!/bin/bash")
    println(f, "")
    #println(f, "#\$ -M mhe@nd.edu")
    #println(f, "#\$ -m abe") # Send mail when job begins, ends and aborts
    #println(f, "#\$ -pe smp 12") # Specify parallel environment and legal core size ("smp" for shared memory)
    println(f, "#\$ -q long") # Specifiy queue ("long" has 14-day run-time limit)
    println(f, "#\$ -N optimize_model") # Specify job name
    println(f, "#\$ -t 1-$n") # Specify number of tasks in array (for job arrays)
    println(f, "")
    println(f, "module load Julia/1.8.5")
end

"""
Generates a job script for running "optimize.jl".
"""
function generate_job_array_optimize(n::Int64=1)
    f_name = "optimize_job.sh"
    f = open(f_name, "w")
    write_job_array_settings_crc(f; n=n)
    #println(f, "julia --project=@. optimize.jl")
    println(f, "julia --project=@. optimize.jl \$SGE_TASK_ID") # for job arrays
    close(f)

    return f_name
end

"""
Generates and submits a job script to run "optimize.jl".
"""
function submit_job_array_optimize(n::Int64=1)
    f_name = generate_job_array_optimize(n)
    run(`qsub $f_name`)
    println("Job array ", f_name, " submitted.")
end





##### NOTE: must run this script from the same directory in which you want to submit the jobs from!

n_runs = 50 # total number of runs/tasks in job array to submit

#submit_job_array_optimize(n_runs)
