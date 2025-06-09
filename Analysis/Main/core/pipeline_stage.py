import logging

class PipelineStage:
    """
    Base class for all pipeline stages.
    """
    def __init__(self, data_manager, config):
        self.data_manager = data_manager
        self.config = config
        self.stage_name = self.__class__.__name__

    def execute(self):
        logging.info(f"--- Executing Stage: {self.stage_name} ---")
        try:
            self._run()
            logging.info(f"--- Stage {self.stage_name} Completed ---")
            return True
        except Exception as e:
            logging.error(f"Error in stage {self.stage_name}: {e}", exc_info=True)
            return False

    def _run(self):
        raise NotImplementedError("Each stage must implement the _run method.")

    def _load_required_data(self, required_data_keys_paths_loaders):
        """
        Helper to load multiple data items if not already in DataManager.
        :param required_data_keys_paths_loaders: List of tuples:
               (data_key, file_path_attr_in_config, load_method_name_in_datamanager, is_object=False)
        :return: True if all loaded successfully or already present, False otherwise.
        """
        all_loaded_successfully = True
        for key, path_attr, load_method_name, is_object in required_data_keys_paths_loaders:
            data_check = self.data_manager.get_object(key) if is_object else self.data_manager.get_data(key)
            if data_check is not None:
                logging.debug(f"Required data/object '{key}' for {self.stage_name} already in memory.")
                continue  # Move to the next required item
            else:
                logging.warning(f"Required data/object '{key}' for {self.stage_name} not found in memory. Attempting to load from file.")
                try:
                    file_path = getattr(self.config, path_attr)
                    load_method = getattr(self.data_manager, load_method_name)
                except AttributeError as e:
                    logging.error(f"Configuration or DataManager attribute error for '{key}' in {self.stage_name}: {e}. Cannot load.")
                    all_loaded_successfully = False
                    continue # Skip to the next item, this one can't be loaded

                if not file_path:
                    logging.error(f"File path for '{key}' resolved to an empty string or None from config attribute '{path_attr}' in {self.stage_name}.")
                    all_loaded_successfully = False
                    continue

                if not load_method(key, file_path):
                    # The loader method itself should log the specifics of the file loading failure.
                    logging.error(f"Failed to load required data/object '{key}' for {self.stage_name} from {file_path}. Stage may not run correctly.")
                    all_loaded_successfully = False
                # Re-fetch after attempting load
                else:
                    data_check = self.data_manager.get_object(key) if is_object else self.data_manager.get_data(key)
                    if data_check is None:
                        logging.error(f"CRITICAL: DataManager method '{load_method_name}' reported success for '{key}' but item is still None in {self.stage_name}.")
                        all_loaded_successfully = False
                    else:
                        logging.info(f"Successfully loaded '{key}' for {self.stage_name} from file.")
        return all_loaded_successfully