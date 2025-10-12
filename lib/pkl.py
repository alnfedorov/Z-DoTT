import pickle
from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class PklData[T]:
    path: Path = field(default_factory=Path)

    def load(self) -> T:
        with open(self.path, 'rb') as stream:
            return pickle.load(stream)

    def dump(self, data: T):
        self.path.parent.mkdir(parents=True, exist_ok=True)

        # First - dump to a temporary file to avoid corrupting existing data
        tmp_path = self.path.with_suffix(".pkl.tmp")
        with open(tmp_path, 'wb') as stream:
            pickle.dump(data, stream, protocol=pickle.HIGHEST_PROTOCOL)

        # Then - rename the temporary file to the original file
        if self.path.exists():
            self.path.unlink()
        tmp_path.rename(self.path)
