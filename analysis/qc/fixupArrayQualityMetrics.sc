import ammonite.ops._
import scala.util.control.NonFatal
val fileToFind = "arrayQualityMetrics.js"
@main
def mainz(wd: Path = pwd) {
    // val matchingFiles = ls.rec(wd).filter(
    //   f => grep(fileToFind, f) match {
    //     case Some(_) => true
    //     case None => false
    //   }
    // )
    // Much faster, Ammonite equivalent, to the snippet above
    val matchingFiles = ls.rec.iter! wd |? grep! fileToFind
    val patches = matchingFiles.map(f => getContentsOfPatchedFileIfBroken((read.lines! f).zipWithIndex))
    val toPatch = matchingFiles.zip(patches.toSeq)

    toPatch.map {
      case (file, Left(msg)) => Left(s"${file}: ${msg}")
      case (file, Right(fileContents)) => backupAndWriteFile(file, fileContents)
    } foreach {
      case Left(msg) => println(msg)
      case Right(msg) => println(msg)
    }
}
def backupAndWriteFile(file: Path, fileContents: Vector[String]): Either[String, String] = {
  try {
    mv(file, file/up/s"${fileToFind}.bak")
    write(file, fileContents.mkString("\n"))
    Right(s"${file}: file backed up and patched successfully")
  } catch {
    case NonFatal(_) => Left(s"${file}: something went wrong!")
  }
}
def getContentsOfPatchedFileIfBroken(fileContents: Vector[(String, Int)]): Either[String, Vector[String]] = {
  BrokenIfMatcher.findBrokenIfRange(fileContents) match {
    case None => Left("Nothing to do")
    case Some((x, y)) => Right(patchFile(fileContents, x, y))
  }
}
def patchFile(fileContents: Vector[(String, Int)], brokenCodeStart: Int, brokenCodeEnd: Int): Vector[String] = {
  (fileContents.take(brokenCodeStart).map(_._1)
    :+ "            if (!selector) {"
    :+ "                i++;"
    :+ "                continue;"
    :+ "            }") ++ (fileContents.takeRight(fileContents.size - brokenCodeEnd - 1).map(_._1))
}
object BrokenIfMatcher {
  val ifLine = """^\s*if\s*\(!selector\)\s*$"""
  val thenLine1 = """^\s*i\+\+;\s*$"""
  val thenLine2 = """^\s*continue;.*$"""
  val ifThenLine = """^\s*if\s*\(!selector\)\s*continue;.*$"""

  def findBrokenIfRange(fileContents: Vector[(String, Int)]): Option[(Int, Int)] = {
    val threeLineWindowIt = fileContents.sliding(3, 1).map { case Vector(e1, e2, e3) => (e1, e2, e3) }
    val threeMatchIt = threeLineWindowIt.filter(l3 => l3._1._1.matches(ifLine) && l3._2._1.matches(thenLine1) && l3._3._1.matches(thenLine2))
    if (threeMatchIt.hasNext) {
      val matchedCode = threeMatchIt.next
      return Some(matchedCode._1._2, matchedCode._3._2)
    }
    val twoLineWindowIt = fileContents.sliding(2, 1).map { case Vector(e1, e2) => (e1, e2) }
    val twoMatchIt = twoLineWindowIt.filter(l2 => l2._1._1.matches(ifLine) && l2._2._1.matches(thenLine2))
    if (twoMatchIt.hasNext) {
      val matchedCode = twoMatchIt.next
      return Some(matchedCode._1._2, matchedCode._2._2)
    }
    val oneLinerMatch = fileContents.filter(l => l._1.matches(ifThenLine))
    if (oneLinerMatch.size == 1) return Some(oneLinerMatch.head._2, oneLinerMatch.head._2)
    return None
  }
}
