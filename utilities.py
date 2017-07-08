import time
from functools import wraps

def fn_timer(function):
    @wraps(function)
    def function_timer(*args, **kwargs):
        t0 = time.time()
        result = function(*args, **kwargs)
        t1 = time.time()
        totaltime = t1 - t0
        if totaltime < 60:
            print ("Total time running %s: %s seconds" %
                (function.func_name, str(totaltime))
                )
        elif totaltime < 3600:
            print ("Total time running %s: %s minutes" %
               (function.func_name, str(totaltime / 60.))
               )
        elif totaltime < 3600 * 24:
            print ("Total time running %s: %s hours and %s minutes" %
               (function.func_name, str(totaltime / 3600),str((totaltime % 3600) / 60.))
               )
        else:
            print ("Total time running %s: %s days,%s hours, and %s minutes" %
               (function.func_name, str(totaltime / (3600 * 24)),str((totaltime % (3600 * 24)) / 3600),str((totaltime % 3600) / 60.))
               )
        return result
    return function_timer