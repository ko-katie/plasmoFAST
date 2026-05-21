import { loadReference } from './reference';
import { getTextStream } from './parser';
import { getCategory } from './classify';
import type { AnalysisResult, StrainResult } from './types';

const DEFAULT_REF_URL = new URL('../reference/25mer_rc_list.tsv', import.meta.url);

self.onmessage = async (event: MessageEvent<{ file: File; referenceUrl?: string }>) => {
  const { file, referenceUrl } = event.data;

  try {
    const ref = await loadReference(referenceUrl ?? DEFAULT_REF_URL.href);

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
          // sequence line — slide 25-mer window
          for (let i = 0; i <= line.length - 25; i++) {
            const kmer = line.slice(i, i + 25);
            const specPos = ref.specKmers.get(kmer);
            if (specPos !== undefined) {
              ref.positions.get(specPos)!.specCount++;
            } else {
              const nonspecPos = ref.nonspecKmers.get(kmer);
              if (nonspecPos !== undefined) {
                ref.positions.get(nonspecPos)!.nonspecCount++;
              }
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
