// Rasterises one or more generated SVG diagrams so that layout problems (clipped labels,
// overlap, overflow) can be inspected without a browser.
// Usage: node scripts/preview-diagram.mjs content/assets/diagrams/<group>/<lang>/<name>.svg ...
import fs from "node:fs"
import path from "node:path"
import { fileURLToPath } from "node:url"
import sharp from "sharp"

const scriptDir = path.dirname(fileURLToPath(import.meta.url))
const repoRoot = path.resolve(scriptDir, "..")
const outputDir = path.join(repoRoot, ".diagram-preview")
fs.mkdirSync(outputDir, { recursive: true })

const inputs = process.argv.slice(2)
if (inputs.length === 0) {
  console.error("Pass at least one SVG path.")
  process.exit(1)
}

for (const input of inputs) {
  const absolute = path.resolve(repoRoot, input)
  const parts = path.relative(path.join(repoRoot, "content", "assets", "diagrams"), absolute)
  const name = parts
    .split(path.sep)
    .join("-")
    .replace(/\.svg$/, ".png")
  const target = path.join(outputDir, name)
  await sharp(absolute, { density: 132 }).png().toFile(target)
  console.log(target)
}
