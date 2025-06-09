import logging

class Pipeline:
    """
    Manages and executes a sequence of PipelineStage objects.
    """
    def __init__(self, data_manager, config):
        self.data_manager = data_manager
        self.config = config
        self.stages = []

    def add_stage(self, stage_class_or_instance):
        if isinstance(stage_class_or_instance, type): # if it's a class
            stage_instance = stage_class_or_instance(self.data_manager, self.config)
        else: # if it's already an instance
            stage_instance = stage_class_or_instance
            
        self.stages.append(stage_instance)
        logging.info(f"Added stage: {stage_instance.stage_name}")

    def run(self):
        logging.info("=== Starting Data Analysis Pipeline ===")
        for stage in self.stages:
            if not stage.execute():
                logging.error(f"Pipeline execution halted due to error in stage: {stage.stage_name}")
                break
        logging.info("=== Data Analysis Pipeline Complete ===")