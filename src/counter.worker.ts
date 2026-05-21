self.addEventListener('error', (e) => {
  self.postMessage({ type: 'error', message: `Worker load error: ${e.message} (${e.filename}:${e.lineno})` });
});

import { parseReference, loadReference } from './reference';
import { getTextStream } from './parser';
import { getCategory } from './classify';
import type { AnalysisResult, StrainResult } from './types';
import defaultRefText from '../reference/25mer_rc_list.tsv';

self.onmessage = async (event: MessageEvent<{ file: File; referenceUrl?: string }>) => {
  const { file, referenceUrl } = event.data;

  try {
    const ref = referenceUrl
      ? await loadReference(referenceUrl)
      : parseReference(defaultRefText);

    const reader = getTextStream(file).getReader();
    let leftover = '';
    let lineIndex = 0;
    let readCount = 0;
    let bytesRead = 0;
    const totalBytes = file.size;

    while (true) {
      const { done, value } = await reader.read();

      const text = leftover + (value ?? '');
      const lines = text.split('\n');
      leftover = lines.pop() ?? '';

      for (const line of lines) {
        if (line.trim() === '') continue;
        if (lineIndex % 4 === 1) {
          for (let i = 0; i <= line.length - 25; i++) {
            const kmer = line.slice(i, i + 25);
            const hit = ref.kmers.get(kmer);
            if (hit !== undefined) {
              if (hit.spec) hit.entry.specCount++;
              else hit.entry.nonspecCount++;
            }
          }
          readCount++;
          if (readCount % 10_000 === 0) {
            if (value) bytesRead += value.length;
            self.postMessage({ type: 'progress', bytesRead, totalBytes });
          }
        }
        lineIndex++;
      }

      if (done) break;
    }

    const result: AnalysisResult = {};
    for (const strain of ref.strains) {
      result[strain] = { specific: 0, nonspecific: 0, mixed: 0, lowCoverage: 0 };
    }
    for (const { specCount, nonspecCount, strain } of ref.positions.values()) {
      (result[strain] as StrainResult)[getCategory(specCount, nonspecCount)]++;
    }

    self.postMessage({ type: 'result', data: result });
  } catch (err) {
    self.postMessage({ type: 'error', message: String(err) });
  }
};
