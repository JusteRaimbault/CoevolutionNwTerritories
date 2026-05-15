package org.nlogo.extensions.transferentropy

import java.io

import org.nlogo.nvm.ExtensionContext

import scala.jdk.CollectionConverters.IteratorHasAsScala
import scala.reflect.Selectable.reflectiveSelectable
import org.nlogo.api.ScalaConversions.RichAny

import org.nlogo.core.LogoList
import org.nlogo.api._
import org.nlogo.core.{ NumberParser, Syntax }
import org.nlogo.api.Reporter
import org.nlogo.api.Command
import org.nlogo.core.Syntax._

import infodynamics.measures.continuous.kraskov.TransferEntropyCalculatorKraskov
import infodynamics.utils.MatrixUtils
import infodynamics.utils.EmpiricalMeasurementDistribution

class TransferEntropyExtension extends DefaultClassManager {

  class TransferEntropy extends Reporter {
    override def getSyntax = Syntax.reporterSyntax(right = List(ListType, ListType, NumberType, NumberType), ret = NumberType)
    def report(args: Array[Argument], context: Context): AnyRef = {
      val X: Array[Double] = try args(0).getList.toArray.map(_.asInstanceOf[Double])
      catch {case e: LogoException => throw new ExtensionException(e.getMessage)}
      val Y: Array[Double] = try args(1).getList.toArray.map(_.asInstanceOf[Double])
      catch {case e: LogoException => throw new ExtensionException(e.getMessage)}

      val k = args(2).getIntValue
      val tau = args(3).getIntValue

      val teCalc: TransferEntropyCalculatorKraskov = new TransferEntropyCalculatorKraskov
      teCalc.setProperty("k", k.toString)
      teCalc.setProperty("delay", tau.toString)
      teCalc.initialise(1)
      teCalc.setObservations(X, Y)
      val res = teCalc.computeAverageLocalOfObservations()
      res.toLogoObject
    }
  }

  class TransferEntropyEnsemble extends Reporter {
    override def getSyntax = Syntax.reporterSyntax(right = List(ListType, ListType, NumberType, NumberType, NumberType), ret = ListType)
    def report(args: Array[Argument], context: Context): AnyRef = {
      val X: Array[Array[Double]] = try args(0).getList.toArray.map(_.asInstanceOf[LogoList].toArray.map(_.asInstanceOf[Double]))
      catch {case e: LogoException => throw new ExtensionException(e.getMessage)}
      val Y: Array[Array[Double]] = try args(1).getList.toArray.map(_.asInstanceOf[LogoList].toArray.map(_.asInstanceOf[Double]))
      catch {case e: LogoException => throw new ExtensionException(e.getMessage)}

      val k = args(2).getIntValue
      val tau = args(3).getIntValue
      val B = args(4).getIntValue

      val teCalc: TransferEntropyCalculatorKraskov = new TransferEntropyCalculatorKraskov
      teCalc.setProperty("k", k.toString)
      teCalc.setProperty("delay", tau.toString)
      teCalc.initialise(1)
      teCalc.startAddObservations()
      X.zip(Y).foreach{case (x,y) =>
        teCalc.addObservations(x, y)
      }
      teCalc.finaliseAddObservations()
      val rawTE: Double = teCalc.computeAverageLocalOfObservations()
      val surrogateDist: EmpiricalMeasurementDistribution = teCalc.computeSignificance(B)
      val ensembleBias: Double = surrogateDist.getMeanOfDistribution
      val TE: Double = rawTE - ensembleBias
      val tScore: Double = surrogateDist.getTSscore
      List(TE, rawTE, tScore).toLogoObject
    }
  }

  override def load(manager: PrimitiveManager) = {
    manager.addPrimitive("transfer-entropy", new TransferEntropy)
    manager.addPrimitive("transfer-entropy-ensemble", new TransferEntropyEnsemble)
  }



}
