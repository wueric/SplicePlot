def check_availability(module_list):
    module_count = 0
    for module_name in module_list:
        try:
            __import__(module_name)
            module_count += 1
        except ImportError:
            print '{0} is either not installed or not in the $PYTHONPATH'.format(module_name)
            pass

    if module_count == len(module_list):
        print 'All modules available'
    return


if __name__ == '__main__':
    module_list = ['numpy','scipy','matplotlib','pandas','pysam']
    check_availability(module_list)    
