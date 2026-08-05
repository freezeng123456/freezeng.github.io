// Prints a title/venue comparison between scripts/data/tao-zhou-publications.json and
// the Crossref record for each recorded DOI, so that a wrong match is easy to spot.
import fs from "node:fs"
import path from "node:path"
import { fileURLToPath } from "node:url"

const scriptDir = path.dirname(fileURLToPath(import.meta.url))
const dataPath = path.join(scriptDir, "data", "tao-zhou-publications.json")
const data = JSON.parse(fs.readFileSync(dataPath, "utf8"))

const normalize = (value) =>
  String(value)
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, " ")
    .trim()

const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms))

for (const paper of data.papers) {
  if (!paper.doi) continue
  const response = await fetch(`https://api.crossref.org/works/${encodeURIComponent(paper.doi)}`, {
    headers: { "User-Agent": "freezeng-kb (verification)" },
  })
  if (!response.ok) {
    console.log(`${paper.n}\tHTTP ${response.status}\t${paper.doi}`)
    await sleep(1200)
    continue
  }
  const item = (await response.json()).message
  const remoteTitle = (Array.isArray(item.title) ? item.title[0] : item.title) ?? ""
  const remoteVenue = (item["container-title"] ?? [])[0] ?? ""
  const remoteYear = item.issued?.["date-parts"]?.[0]?.[0] ?? "?"
  const same = normalize(remoteTitle) === normalize(paper.title)
  console.log(
    `${paper.n}\t${same ? "OK  " : "DIFF"}\t${paper.doi}\t${remoteYear}\t${remoteVenue}\n\tlocal : ${paper.title}\n\tremote: ${remoteTitle}`,
  )
  await sleep(1200)
}
