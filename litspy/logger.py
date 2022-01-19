import logging
import pathlib


class Logger:
    """
    class for creating and initialising the logger
    """
    def __init__(self):
        """
        initialise the dictionary of logging levels
        """
        self.log_levels_dict = {10: 'DEBUG', 20: 'INFO', 30: 'WARNING', 40: 'ERROR', 50: 'CRITICAL'}

    @staticmethod
    def initialise_logging_to_file(logger, datestamp):
        """
        initialise logging to a file with a file name (datestamp) and log formatting

        :param logging.Logger logger: logger object
        :param str datestamp: string representing date and time the process started
        :return: logger object with configured file handler
        :rtype: logging.Logger
        """
        # name the file with the datestamp, or give a default name if no datestamp
        if datestamp:
            file_path = pathlib.Path(f'{datestamp}.log')
        else:
            file_path = pathlib.Path('litspy.log')

        # create the file handler
        file_h = logging.FileHandler(file_path, encoding='utf-8')

        # create and apply the formatter to the file handler
        file_fmt = logging.Formatter('%(asctime)s \t%(levelname)s:\t%(message)s', datefmt='%Y/%m/%d %H:%M:%S')
        file_h.setFormatter(file_fmt)

        # set file logging level to info
        file_h.setLevel(logging.INFO)

        # add the file handler to the logger object
        logger.addHandler(file_h)
        return logger

    @staticmethod
    def initialise_logging_to_console(int_level, logger):
        """
        initialise logging at the console level

        :param int int_level: logging level (int)
        :param logging.Logger logger: logger object
        :return: logger object with configured console handler
        :rtype: logging.Logger
        """
        # initialise console handler and set the console logging level
        cons_h = logging.StreamHandler()
        cons_h.setLevel(int_level)

        # set console formatter and add it to the console handler
        cons_fmt = logging.Formatter('%(asctime)s|%(levelname)s: %(message)s', datefmt='%H:%M:%S')
        cons_h.setFormatter(cons_fmt)

        # add the handler to the logger and return the logger
        logger.addHandler(cons_h)
        return logger

    def get_logging_integer(self, supplied_level):
        """
        determine the logging level by comparing the supplied logging level to the dictionary of log levels, return the
        logging level and whether it is the default
         
        :param str or int supplied_level: word or number describing logging level
        :return: supplied or default logging level, bool for whether the default level was used (True = default used)
        :rtype: tuple[int, bool]
        """
        # initialise
        default_used = True
        level = 30

        # convert the logging level to an integer if possible
        try:
            supplied_level = int(supplied_level)
        except ValueError:
            pass

        # if the supplied level is a valid logging level value, return it in the accepted integer format
        if isinstance(supplied_level, int):
            if supplied_level in self.log_levels_dict.keys():
                default_used = False
                level = supplied_level
            elif supplied_level*10 in self.log_levels_dict.keys():
                default_used = False
                level = supplied_level*10
        elif supplied_level.upper() in self.log_levels_dict.values():
            int_level = getattr(logging, supplied_level.upper())
            default_used = False
            level = int_level

        return level, default_used

    def get_logging_level(self, supplied_level):
        """
        identify and return the supplied logging level, or the default level if an unexpected value or no value was
        provided

        :param str or int or list supplied_level: logging level provided
        :return: logging level (int), default_used (bool)
        :rtype: tuple[int, bool]
        """
        # initialise to 'warning'
        level = 30
        default_used = True

        # if a warning level was specified
        if supplied_level:
            # if the specified level is a single-element list, then extract the level from the list
            if isinstance(supplied_level, list):
                if len(supplied_level) == 1:
                    level_from_list = supplied_level[0]
                    # convert the supplied value to an integer logging level, and return whether the default was used
                    level, default_used = self.get_logging_integer(level_from_list)
            else:
                # convert the supplied value to an integer logging level, and return whether the default was used
                level, default_used = self.get_logging_integer(supplied_level)
        return level, default_used

    def initialise_logger(self, level=None, logfile=None, datestamp=None):
        """
        initialise the logger to be used throughout the package

        :param str or int or list level: supplied console logging level
        :param bool logfile: whether to log to file
        :param str datestamp: string representing date and time the process started
        :return: logger object
        :rtype: logging.Logger
        """
        # initialise the logger
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)

        # if logfile has been specified, initialise logging to file
        if logfile:
            logger = self.initialise_logging_to_file(logger, datestamp)

        # get and set the console logging level
        level, default_used = self.get_logging_level(level)
        logger = self.initialise_logging_to_console(level, logger)

        # log relevant messages about logging initialisation
        logger.info(f"Console logging level set to '{self.log_levels_dict[level]}'")

        if not logfile:
            logger.info(f"Logging to console only. To turn on logging to file, run the command again with the -f flag")

        if default_used:
            logger.warning("No valid console logging level supplied; using 'WARNING' level by default")

        return logger
