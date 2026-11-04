"""Small subprocess fixture; no optional native packages or scientific fitting."""
import os
import time
from pathlib import Path
import numpy as np
from codedcal.workflow_validation import Catalog, WorkflowResult


class Fixture:
    def __init__(self, config):
        self.config = config

    def run(self, case, policy):
        assert 'challenge_truth' not in case.metadata
        if policy == 'fail':
            raise ValueError('Deliberate child failure')
        if policy == 'wait':
            time.sleep(60)
        if self.config.get('marker'):
            with Path(self.config['marker']).open('a') as stream:
                stream.write(str(os.getpid()) + '\n')
        return WorkflowResult(Catalog([[0, 2, 3]], [0]), np.zeros((len(case.images), 3)), 1,
            dict(actual_workflow_executions=1, actual_calibration_executions=1, actual_decoder_executions=1,
                 calibration_seconds=.01, decoder_seconds=.02, wall_seconds=.03, exact_result_cache_hit=False))


def make_adapter(config=None):
    return Fixture(config or {})
