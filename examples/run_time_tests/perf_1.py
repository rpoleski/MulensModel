"""
Run benchmarks using perf package. We run the same calculations using
MulensModel, pyLIMA, and numpy.
"""
import perf


def read_and_simplify(file_name):
    """
    Read file and combine it in a single line.
    """
    with open(file_name) as in_file:
        lines = [line for line in in_file.readlines() if line.rstrip()]
        text = ''.join(lines).replace('\n', '; ')
    return text


kwargses = []
kwargses.append(dict(name='MM', setup=read_and_simplify('MM_setup_1.py'), 
    stmt='event.get_chi2()'))
kwargses.append(dict(name='pyLIMA', 
    setup=read_and_simplify('pyLIMA_setup_1.py'), 
    stmt='chi2_telescope(your_event, model_1, pyLIMA_parameters)'))
kwargses.append(dict(name='numpy', 
    setup=read_and_simplify('numpy_setup_1.py'),
    stmt=read_and_simplify('numpy_run_1.py')))

n_processes = 10 # 20 is default value.

# Main part is below.
runner = perf.Runner(processes=n_processes)
for kwargs in kwargses:
    runner.timeit(**kwargs)

