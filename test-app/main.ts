import { analyze } from '../src/index';
import type { AnalysisResult, StrainResult } from '../src/index';

const fileInput = document.getElementById('file-input') as HTMLInputElement;
const progressBar = document.getElementById('progress-bar') as HTMLProgressElement;
const progressLabel = document.getElementById('progress-label') as HTMLSpanElement;
const resultsTable = document.getElementById('results-body') as HTMLTableSectionElement;
const errorBox = document.getElementById('error') as HTMLDivElement;
const status = document.getElementById('status') as HTMLParagraphElement;

fileInput.addEventListener('change', async () => {
  const file = fileInput.files?.[0];
  if (!file) return;

  errorBox.textContent = '';
  resultsTable.innerHTML = '';
  progressBar.value = 0;
  status.textContent = 'Analysing…';

  try {
    const result: AnalysisResult = await analyze(file, {
      onProgress: ({ bytesRead, totalBytes }) => {
        const pct = Math.round((bytesRead / totalBytes) * 100);
        progressBar.value = pct;
        progressLabel.textContent = `${pct}%`;
      },
    });

    for (const [strain, counts] of Object.entries(result) as [string, StrainResult][]) {
      const row = resultsTable.insertRow();
      [strain, counts.specific, counts.nonspecific, counts.mixed, counts.lowCoverage]
        .forEach(val => { row.insertCell().textContent = String(val); });
    }

    progressBar.value = 100;
    progressLabel.textContent = '100%';
    status.textContent = 'Done.';
  } catch (err) {
    errorBox.textContent = String(err);
    status.textContent = '';
  }
});
